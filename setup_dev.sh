#! /bin/bash

# Install all development envs for python code. Assumes mamba is installed

for i in workflow/scripts/{python,rmarkdown}/*; do
    base="$(basename "$i")"
    if [ "$base" != "common" ]; then
        mamba env update -f "$i/env.yml" -n "stratomod-$base" --prune
    fi
done

