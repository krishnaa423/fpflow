#!/bin/bash

if [ -z "$1" ]; then
    echo "Need to provide version as first argument."
fi

twine upload dist/fpflow-$1-py3-none-any.whl dist/fpflow-$1.tar.gz