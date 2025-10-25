#!/bin/bash


. ../../setup.sh

$MPIRUN pw.x -in scf.in &> scf.out

