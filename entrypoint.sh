#!/bin/bash -l

set -e

echo `pwd`

echo "Getting run.py"
wget https://raw.githubusercontent.com/mbackenkoehler/kinodata-docking/main/run.py

python run.py
