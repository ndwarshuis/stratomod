#! /bin/bash

# lint all scripts

# usage: lint.sh [ENV_BASE_PREFIX]
#
# where ENV_BASE_PREFIX is an optional location to store each environment
# (meant to be used by the CI/CD pipeline) and should correspond to what was
# used to create the environments using `setup_dev.sh`

eval "$(${CONDA_EXE} shell.bash hook 2> /dev/null)"
anyfailed=false
if [ -z "$1" ]; then
    prefix="$1"
fi

_mypy () {
    # run mypy without the cache to ensure a clean (albeit slow) lint
    mypy --no-incremental --cache-dir=/dev/null "$1"
}

_conda_activate () {
    # use either a local env or a user-installed env
    base="$1"
    if [ -z "$prefix" ]; then
        conda activate stratomod-"$base"
    else
        conda activate -p "$prefix-$base"
    fi
}

for base in bedtools ebm; do
    echo "Testing scripts for env: $base"
    echo ""
    
    _conda_activate "$base"
    
    flake8 "workflow/scripts/python/$base"
    anyfailed=$([ $? -ne 0 ] || $anyfailed)
    
    _mypy "workflow/scripts/python/$base"
    anyfailed=$([ $? -ne 0 ] || $anyfailed)

    echo ""
done

if $anyfailed; then
   exit 1
fi
