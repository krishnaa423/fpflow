#!/bin/bash


. ../../setup.sh

$MPIRUN pw2bgw.x -pd .true. -in wfn.pp.in &> wfn.pp.out

