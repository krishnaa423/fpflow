#!/bin/bash

cd espresso/02-Wfn
mkdir -p struct.save
cd ../../

cd espresso/03-Wfnq
mkdir -p struct.save
cd ../../

cd espresso/04-Wfn_co
mkdir -p struct.save
cd ../../

cd espresso/05-Wfn_fi
mkdir -p struct.save
cd ../../

cd espresso/06-Wfnq_fi
mkdir -p struct.save
cd ../../

ln -nfs ../../01-Density/struct.save/charge-density.hdf5 espresso/02-Wfn/struct.save/charge-density.hdf5
ln -nfs ../../01-Density/struct.save/spin-polarization.hdf5 espresso/02-Wfn/struct.save/spin-polarization.hdf5

ln -nfs ../../01-Density/struct.save/charge-density.hdf5 espresso/03-Wfnq/struct.save/charge-density.hdf5
ln -nfs ../../01-Density/struct.save/spin-polarization.hdf5 espresso/03-Wfnq/struct.save/spin-polarization.hdf5

ln -nfs ../../01-Density/struct.save/charge-density.hdf5 espresso/04-Wfn_co/struct.save/charge-density.hdf5
ln -nfs ../../01-Density/struct.save/spin-polarization.hdf5 espresso/04-Wfn_co/struct.save/spin-polarization.hdf5

ln -nfs ../../01-Density/struct.save/charge-density.hdf5 espresso/05-Wfn_fi/struct.save/charge-density.hdf5
ln -nfs ../../01-Density/struct.save/spin-polarization.hdf5 espresso/05-Wfn_fi/struct.save/spin-polarization.hdf5

ln -nfs ../../01-Density/struct.save/charge-density.hdf5 espresso/06-Wfnq_fi/struct.save/charge-density.hdf5
ln -nfs ../../01-Density/struct.save/spin-polarization.hdf5 espresso/06-Wfnq_fi/struct.save/spin-polarization.hdf5


ln -nfs ../espresso/02-Wfn/wfn.cplx.h5 11-epsilon/WFN.h5
ln -nfs ../espresso/03-Wfnq/wfn.cplx.h5 11-epsilon/WFNq.h5
ln -nfs ../espresso/04-Wfn_co/wfn.cplx.h5 12-sigma/WFN_inner.h5
ln -nfs ../espresso/04-Wfn_co/rho.real 12-sigma/RHO
ln -nfs ../espresso/04-Wfn_co/vxc.dat 12-sigma/vxc.dat


ln -nfs ../11-epsilon/eps0mat.h5 12-sigma/eps0mat.h5
ln -nfs ../11-epsilon/epsmat.h5 12-sigma/epsmat.h5
ln -nfs ../espresso/04-Wfn_co/wfn.cplx.h5 13-kernel/WFN_co.h5
ln -nfs ../11-epsilon/eps0mat.h5 13-kernel/eps0mat.h5
ln -nfs ../11-epsilon/epsmat.h5 13-kernel/epsmat.h5
ln -nfs ../espresso/04-Wfn_co/wfn.cplx.h5 14-absorption/WFN_co.h5
ln -nfs ../espresso/05-Wfn_fi/wfn.cplx.h5 14-absorption/WFN_fi.h5
ln -nfs ../espresso/06-Wfnq_fi/wfn.cplx.h5 14-absorption/WFNq_fi.h5
ln -nfs ../12-sigma/sigma_hp.log 14-absorption/sigma_hp.log
ln -nfs ../12-sigma/eqp1.dat 14-absorption/eqp_co.dat
ln -nfs ../11-epsilon/eps0mat.h5 14-absorption/eps0mat.h5
ln -nfs ../11-epsilon/epsmat.h5 14-absorption/epsmat.h5
ln -nfs ../13-kernel/bsemat.h5 14-absorption/bsemat.h5
