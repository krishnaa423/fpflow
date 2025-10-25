#!/bin/bash


. ../../setup.sh

cp -f ../01-Density/struct.save/data-file-schema.xml struct.save/data-file-schema.xml
$MPIRUN pw.x -in wfn.in &> wfn.out

