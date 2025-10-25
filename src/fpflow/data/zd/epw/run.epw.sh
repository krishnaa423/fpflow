#!/bin/bash

. ../setup.sh

###Move the KS wave functions for gauge consistency###
rm -rf struct.save 
ln -sf ../espresso/05-Wfn_fi/struct.save struct.save 

###Move all the BSE eigenvectors to the current folder###
if [ -d "eigv" ]; then
  echo -e "BSE eigenvector exists"
else
  echo -e "Start copying BSE eigenvectors"

  mkdir eigv
  for ((i=1;i<=$NQPT;i++))
  do
      j=$((i-1))
      cd eigv
      mkdir q_$j
      cd q_$j
      ln -sf ../../../finiteQ_grid/Q$i/14-absorption/eigenvectors.h5 .
      cd ../../
  done

  echo -e "Done copy BSE eigenvectors"
fi
##########################################################

$MPIRUN_EPW  epw.x -npool $NPOOL  -inp epw1.in > epw1.out
echo -e "Done xctph \n\n"
$MPIRUN_EPW  epw.x -npool $NPOOL  -inp epw2.in > epw2.out
echo -e "Done ste \n\n"
$MPIRUN_EPW  epw.x -npool $NPOOL  -inp epw3.in > epw3.out
echo -e "Done ste eh densities xsf \n\n"
