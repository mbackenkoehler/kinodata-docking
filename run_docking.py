from typing import List
import shutil
import requests as req
import subprocess
import pathlib
import os
import multiprocessing
import time
import psutil
import sys
import signal
import random
from functools import cached_property
from collections import namedtuple

import pandas as pd

HERE = pathlib.Path(".").absolute()


class DockingJob:
    def __init__(self, task):
        self.ident = task.ident
        self.protein_filepath = self.output_dir / "protein.pdb"
        shutil.copy2(task.protein, self.protein_filepath)
        self.smiles = task.smiles
        # maybe the docking was done in a previous run
        if not self.success:
            self._start()
        self.start_time = time.time()

    def _start(self):
        out = open(self.output_dir / "run.log", "w")
        err = open(self.output_dir / "run.err", "w")
        self.process = subprocess.Popen(
            [
                "conda",
                "run",
                "--no-capture-output",
                "-n",
                "kinoml",
                "python",
                "docking.py",
                str(self.ident),
                str(self.protein_filepath),
                str(self.smiles),
                str(self.output_dir),
            ],
            stdout=out,
            stderr=err,
            close_fds=True,
            shell=False,
        )

    @property
    def output_dir(self):
        out_dir = HERE / "docking_pipeline" / "complexes" / str(self.ident)
        out_dir.mkdir(exist_ok=True, parents=True)
        return out_dir

    @property
    def pid(self):
        if hasattr(self, "process"):
            return self.process.pid
        else:
            return -1

    @property
    def success(self):
        """Check if docking results are there for `ident`."""
        output_file = self.output_dir / "docking.csv"
        ligand_file = self.output_dir / f"{self.ident}_ligand.pdb"
        success = output_file.exists() and ligand_file.exists()
        return success

    @property
    def running(self):
        return self.success or self.process.poll()

    @property
    def memory(self):
        """overall process tree's memory [gb]"""
        try:
            return sum(
                child.memory_info().rss / 1024**3
                for child in psutil.Process(self.pid).children(recursive=True)
            )
        except psutil.NoSuchProcess:
            return 0

    @property
    def runtime(self):
        """runtime in minutes"""
        return (time.time() - self.start_time) / 60

    def suicide(
        self, sig=signal.SIGTERM, include_parent=True, timeout=None, on_terminate=None
    ):
        """Kill a process tree (including grandchildren) with signal
        "sig" and return a (gone, still_alive) tuple.
        "on_terminate", if specified, is a callback function which is
        called as soon as a child terminates.
        """
        assert self.pid != os.getpid(), "I won't kill the master"
        try:
            os.remove(self.protein_filepath)
        except FileNotFoundError:
            pass
        if self.pid < 0:
            return
        try:
            parent = psutil.Process(self.pid)
        except psutil.NoSuchProcess:
            return
        children = parent.children(recursive=True)
        if include_parent:
            children.append(parent)
        for p in children:
            try:
                p.send_signal(sig)
            except psutil.NoSuchProcess:
                pass
        gone, alive = psutil.wait_procs(
            children, timeout=timeout, callback=on_terminate
        )
        return (gone, alive)


class DockingScheduler:
    def __init__(
        self,
        tasks,
        capacity=os.cpu_count(),
        proc_mem_limit=10,
        timeout=10,
        total_mem_start_limit=50,
        output_dir=HERE / "docking_pipeline",
    ):
        """
        Parameters
        ----------
        tasks: List[Task]
            docking tasks
        capacity: int
            the maximum number of concurrent docking jobs
        proc_mem_limit: int
            memory limit per job in gb
        timeout: int
            job timeout in minutes
        total_mem_start_limit: int
            limit in percentage on memory above which no new jobs are started
        """
        self.capacity = capacity
        self.proc_mem_limit = proc_mem_limit
        self.timeout = timeout
        self.total_mem_start_limit = total_mem_start_limit
        self.running = list()
        self.waitlist = tasks
        self.output_dir = output_dir

    @property
    def failure_file(self):
        return self.output_dir / "failures.csv"

    def print_status(self):
        print(
            f"\r|waitlist| = {len(self.waitlist)} |running| = {len(self.running)}  ",
            flush=True,
            end="\r",
        )

    def run(self):
        while len(self.waitlist) > 0 or len(self.running) > 0:
            self.print_status()
            time.sleep(1)  # busy wait...
            self.cleanup_running()

            # check overall memory usage
            if psutil.virtual_memory()[2] > self.total_mem_start_limit:
                continue

            self.start_dockings()

    def start_dockings(self):
        for _ in range(self.capacity - len(self.running)):
            if len(self.waitlist) == 0:
                return
            task = self.waitlist.pop()
            self.running.append(DockingJob(task))

    def cleanup_running(self):
        # clean self.running processes
        still_running = list()
        for i, job in enumerate(self.running):
            # check for completion
            if job.running:
                print(f"{job.ident} is done")
                job.suicide()  # make sure it's dead
                continue

            # check for timeout and out-of-memory
            if job.memory > self.proc_mem_limit:
                with open(self.failure_file, "a") as f:
                    f.write(f"{job.ident},memory\n")
                job.suicide()
                continue
            if job.runtime > self.timeout:
                with open(self.failure_file, "a") as f:
                    f.write(f"{job.ident},timeout\n")
                job.suicide()
                continue
            still_running.append(i)
        self.running = [job for i, job in enumerate(self.running) if i in still_running]


Task = namedtuple("Task", "ident protein smiles")


class TemplateData:
    def __init__(
        self,
        kinodata_path="activities-chembl31.csv",
        similar_pdb_path="docking_pipeline/most_similar.csv",
    ):
        self.kinodata_path = kinodata_path
        self.similar_pdb_path = similar_pdb_path

    @cached_property
    def kinodata(self):
        # activities.activity_id,assays.chembl_id,target_dictionary.chembl_id,molecule_dictionary.chembl_id,molecule_dictionary.max_phase,activities.standard_type,activities.standard_value,activities.standard_units,compound_structures.canonical_smiles,compound_structures.standard_inchi,component_sequences.sequence,assays.confidence_score,docs.chembl_id,docs.year,docs.authors,UniprotID
        return pd.read_csv(self.kinodata_path, index_col="activities.activity_id")

    @cached_property
    def similar_pdbs(self):
        # activities.activity_id,similar.ligand_pdb,similar.complex_pdb,similar.chain
        return pd.read_csv(
            self.similar_pdb_path, index_col="activities.activity_id"
        ).head(20)


def prepare_tasks(
    data: TemplateData, output_dir=HERE / "docking_pipeline"
) -> List[Task]:
    klifs_structures_file = output_dir / "klifs_structures.csv"
    if klifs_structures_file.exists():
        structures = pd.read_csv(
            klifs_structures_file, index_col="activities.activity_id"
        )
    else:
        structures = None
    print("-> populate waitlist")
    tasks, idents, ids = list(), list(), list()
    for ident, row in data.similar_pdbs.iterrows():
        if structures is None or ident not in structures.index:
            try:
                structure_ID = get_klifs_structure_id(
                    row["similar.complex_pdb"],
                    row["similar.chain"],
                    row["similar.ligand_pdb"],
                )
            except ValueError:
                continue
            idents.append(ident)
            ids.append(structure_ID)
        else:
            structure_ID = structures.loc[ident, "similar.klifs_structure_id"]
        protein_file = get_klifs_protein(
            structure_ID, output_path=output_dir / "KLIFS_proteins"
        )
        task = Task(
            ident,
            protein_file,
            data.kinodata.loc[ident, "compound_structures.canonical_smiles"],
        )
        tasks.append(task)
    structure_ids = pd.DataFrame(
        {"activities.activity_id": idents, "similar.klifs_structure_id": ids}
    ).set_index("activities.activity_id")
    if structures is not None:
        structure_ids = pd.concat([structure_ids, structures])

    structure_ids.to_csv(klifs_structures_file)

    return tasks


def get_klifs_structure_id(pdb_id: str, chain: str, ligand_id: str):
    """
    Get the complex PDB from KLIFS.

    Parameters
    ----------
    pdb_id: str
        PDB ID.
    chain: str
        the chain.
    ligand_id: str
        ligand PDB ID

    Returns
    -------
    structure_ID: int
        KLIFS structure id for the most similar
    """
    resp = req.get(
        "https://klifs.net/api_v2/structures_pdb_list", {"pdb-codes": pdb_id}
    )
    klifs_info = None
    resp.raise_for_status()
    for info in resp.json():
        if (
            str(info["chain"]).upper() == str(chain).upper()
            and str(info["ligand"]).upper() == str(ligand_id).upper()
        ):
            if klifs_info is None:
                klifs_info = info
            elif klifs_info["quality_score"] < info["quality_score"]:
                klifs_info = info
    if klifs_info is None:
        raise ValueError(f"not found pdb:{pdb_id} chain:{chain} lig_pdb:{ligand_id}")
    return klifs_info["structure_ID"]


def get_klifs_protein(structure_ID: int, output_path=HERE):
    """
    Get the complex mol2 from KLIFS.

    Parameters
    ----------
    structure_ID: int
        KLIFS structure ID
    path: Path, optional
        folder to store the structure in.

    Returns
    -------
    file_path: Path
        path of the mol2 file.
    """
    filename = output_path / f"{structure_ID}.pdb"
    if filename.exists():
        return filename
    pathlib.Path(output_path).mkdir(exist_ok=True, parents=True)
    resp = req.get(
        "https://klifs.net/api_v2/structure_get_pdb_complex",
        {"structure_ID": structure_ID},
    )

    with open(filename, "w") as f:
        f.write(resp.text)

    return filename


if __name__ == "__main__":

    print("-> read data")
    data = TemplateData()
    print("-> prepare docking tasks")
    tasks = prepare_tasks(data)

    scheduler = DockingScheduler(
        tasks,
        capacity=os.cpu_count(),
        proc_mem_limit=1,
        timeout=5,
        total_mem_start_limit=50,
    )

    print("-> start docking")
    scheduler.run()
