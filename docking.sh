#!/bin/bash -l

set -e

# echo "Getting run.py"
# wget https://raw.githubusercontent.com/mbackenkoehler/kinodata-docking/main/run.py
cd /home/michael.backenkoehler/docking

conda run --no-capture-output -n kinoml python docking.py
