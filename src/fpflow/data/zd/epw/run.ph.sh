#!/bin/bash

. ../setup.sh

cp ../pseudos/*.upf .

$MPIRUN  pw.x -inp scf.in > scf.out
echo -e "Done epw scf \n\n"
$MPIRUN  ph.x -inp ph.in > ph.out
echo -e "Done epw ph \n\n"
