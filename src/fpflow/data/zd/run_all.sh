#!/bin/bash

# Run QE and BSE. 
./run.espresso.sh

# Run BSEq. 
cd ./finiteQ_grid
python ./generate_ex_band.py
./run.bse.sh
cd ../

# Run epw. 
cd ./epw
./run.ph.sh
python ./create_save.py
./run.epw.sh
cd ../