#! /bin/bash

# Install all development envs for python code. Assumes mamba is installed

# usage: setup_dev.sh [ENV_BASE_PREFIX]
#
# where ENV_BASE_PREFIX is an optional location to store each environment
# (meant to be used by the CI/CD pipeline)

for i in workflow/scripts/{python,rmarkdown}/*; do
    base="$(basename "$i")"
    if [ "$base" != "common" ]; then
        if [ -z "$1" ]; then
            location=("-n" "stratomod-$base")
        else
            location=("-p" "$1-$base")
        fi
        mamba env update -f "$i/env.yml" "${location[@]}" --prune
        
    fi
done

