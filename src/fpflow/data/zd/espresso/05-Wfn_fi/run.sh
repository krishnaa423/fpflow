#!/bin/bash


. ../../setup.sh

bash wfn.run.sh
bash run-pw2bgw.sh
wfn2hdf.x BIN wfn.cplx wfn.cplx.h5