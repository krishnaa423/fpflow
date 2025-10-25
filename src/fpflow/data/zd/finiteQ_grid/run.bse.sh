#!/bin/bash

. ../setup.sh

for ((i=1; i<=$NQPT; i++))
do
    echo "Starting BSE at the ${i}-th(/$NQPT) Q point..."

    cd Q$i
    bash link.sh

    cd 13-kernel
    bash run.sh
    cd ..

    cd 14-absorption
    bash run.sh
    cd ../13-kernel
    rm bsemat.h5
    cd ../../
done

