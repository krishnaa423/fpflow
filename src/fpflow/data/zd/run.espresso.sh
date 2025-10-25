#!/bin/bash
bash espresso/link.sh
cd espresso/01-Density
bash run.sh
echo -e "Done Density \n\n"
cd ../..

cd espresso/02-Wfn
mkdir -p struct.save
bash run.sh
echo -e "Done Wfn \n\n"
cd ../..

cd espresso/03-Wfnq
mkdir -p struct.save
bash run.sh
echo -e "Done Wfnq \n\n"
cd ../..

cd espresso/04-Wfn_co
mkdir -p struct.save
bash run.sh
echo -e "Done Wfn_co \n\n"
cd ../..

cd espresso/05-Wfn_fi
mkdir -p struct.save
bash run.sh
echo -e "Done Wfn_fi \n\n"
cd ../..

cd espresso/06-Wfnq_fi
mkdir -p struct.save
bash run.sh
echo -e "Done Wfnq_fi \n\n"
cd ../..

cd 11-epsilon
bash run.sh
echo -e "Done epsilon \n\n"
cd ..

cd 12-sigma
bash run.sh
echo -e "Done sigma \n\n"
cd ..

cd 13-kernel
bash run.sh
echo -e "Done kernel \n\n"
cd ..

cd 14-absorption
bash run.sh
echo -e "Done absorption \n\n"
cd ..


