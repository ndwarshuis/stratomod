#! /bin/bash

# lint all scripts

# usage: lint.sh [ENV_BASE_PREFIX]
#
# where ENV_BASE_PREFIX is an optional location to store each environment
# (meant to be used by the CI/CD pipeline) and should correspond to what was
# used to create the environments using `setup_dev.sh`

anyfail=0 # global var (ew...) to track if we have any lint failures

_run_test () {
    "$@"
    local status=$?
    if (( status != 0 )) || (( anyfail == 1 )); then
        anyfail=1
    fi
}

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
        conda activate "./$prefix-$base"
    fi
}

eval "$(${CONDA_EXE} shell.bash hook 2> /dev/null)"

python_root="workflow/scripts/python"

if [ -n "$1" ]; then
    prefix="$1"
fi

for base in bio ebm; do
    echo "Testing scripts for env: $base"
    echo ""
    
    _conda_activate "$base"

    # test the common dir since flake8 doesn't follow imports
    _run_test flake8 "$python_root/common"
    _run_test flake8 "$python_root/$base"
    _run_test _mypy "$python_root/$base"

    echo ""
done

if (( anyfail == 1 )); then
   exit 1
fi
