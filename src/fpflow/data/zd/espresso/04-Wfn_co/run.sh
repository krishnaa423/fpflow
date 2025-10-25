#!/bin/bash


. ../../setup.sh

bash wfn.run.sh
bash run-pw2bgw.sh
parabands.cplx.x &> parabands.out
