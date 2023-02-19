#! /bin/bash

eval "$(${CONDA_EXE} shell.bash hook 2> /dev/null)"
anyfailed=false

_mypy () {
    mypy --no-incremental --cache-dir=/dev/null "$1"
}

echo "Testing bed scripts"
echo ""

conda activate stratomod-bedtools

flake8 workflow/scripts/python/bedtools
anyfailed=$([ $? -ne 0 ] || $anyfailed)

_mypy workflow/scripts/python/bedtools
anyfailed=$([ $? -ne 0 ] || $anyfailed)

echo ""

echo "Testing EBM scripts"
echo ""

conda activate stratomod-ebm

flake8 workflow/scripts/python/ebm
anyfailed=$([ $? -ne 0 ] || $anyfailed)

_mypy workflow/scripts/python/ebm
anyfailed=$([ $? -ne 0 ] || $anyfailed)

if $anyfailed; then
   exit 1
fi
