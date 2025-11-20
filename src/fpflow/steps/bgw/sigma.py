#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.inputs.grammars.namelist import NamelistGrammar
from fpflow.inputs.grammars.bgw import BgwGrammar
import jmespath
from fpflow.io.update import update_dict
from fpflow.io.logging import get_logger
from fpflow.schedulers.scheduler import Scheduler
from importlib.util import find_spec
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.structure.kpts import Kpts
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class BgwSigmaStep(Step):
    @property
    def sigma_update_from_pseudobands(self) -> str:
        return f'''#!/usr/bin/env python

import h5py
import glom 
from fpflow.inputs.inputyaml import InputYaml
from fpflow.steps.bgw.sigma import BgwSigmaStep
from fpflow.structure.qe.qe_struct import QeStruct
import jmespath
import os

inputdict: dict = InputYaml.from_yaml_file('../input.yaml').inputdict
num_pseudobands: int = 0
wfnlink: str = jmespath.search('gw.sigma.wfnlink', inputdict)
parabands_file: str = os.path.join('..', wfnlink, 'wfn_parabands.h5')
with h5py.File(parabands_file, 'r') as f:
    num_pseudobands = f['/parabands/pseudobands/nb_total'][()]
# Qestruct.
max_val_bands: int = int(QeStruct.from_inputdict(inputdict).max_val(
    xc=jmespath.search('scf.xc', inputdict),
    is_soc=jmespath.search('scf.is_spinorbit', inputdict),
))

num_gw_cond_bands: int = num_pseudobands - max_val_bands - 1
glom.assign(inputdict, 'gw.sigma.conv_cond_bands', num_gw_cond_bands)
sig_step = BgwSigmaStep(inputdict=inputdict)
os.chdir('../')
sig_step.create()
os.chdir('./sigma')
'''
    
    @property
    def sigma(self) -> str:
        # Qestruct.
        max_val_bands: int = QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )

        # wfnlink. 
        wfn_link: str = jmespath.search('gw.sigma.wfnlink', self.inputdict)
        wfn_names: list[str] = jmespath.search(f'nscf.list[*].name', self.inputdict)
        wfn_index: int = wfn_names.index(wfn_link)
        self.wfn_options: dict = jmespath.search(f'nscf.list[{wfn_index}]', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('kgrid[0]', self.wfn_options),
                jmespath.search('kgrid[1]', self.wfn_options),
                jmespath.search('kgrid[2]', self.wfn_options),
            ],
            is_reduced=jmespath.search('sym', self.wfn_options),
        )

        sigmadict: dict = {
            'no_symmetries_q_grid': '',
            'number_bands': jmespath.search('gw.sigma.conv_cond_bands', self.inputdict) + max_val_bands,
            'band_index_min': max_val_bands - jmespath.search('gw.sigma.val_bands', self.inputdict) + 1,
            'band_index_max': max_val_bands + jmespath.search('gw.sigma.cond_bands', self.inputdict),
            'degeneracy_check_override': '',
            'screened_coulomb_cutoff': jmespath.search('gw.sigma.ecut', self.inputdict),
            'use_wfn_hdf5': '',
            'kpoints': kpts.sigma_kpts,
        }

        # Update if needed. 
        update_dict(sigmadict, jmespath.search('gw.sigma.args', self.inputdict))

        return BgwGrammar().write(sigmadict)
    
    @property
    def job_sigma(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'gw.sigma.job_info')

        wfnlink_str: str = os.path.join(
            '..',
            jmespath.search('gw.sigma.wfnlink', self.inputdict),
            'wfn_parabands.h5' if jmespath.search('parabands.enabled', self.wfn_options) else 'wfn.h5' 
        )

        unfold_string: str = f'''python -c "from fpflow.analysis.gw.unfold_sigma import unfold_sigma; unfold_sigma()" &> unfold_sigma.out'''

        kgrid_link_folder: str = f"../kgrid_{jmespath.search('kgrid[0]', self.wfn_options)}_{jmespath.search('kgrid[1]', self.wfn_options)}_{jmespath.search('kgrid[2]', self.wfn_options)}_qshift_{0.0}_{0.0}_{0.0}"

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{'python ./script_update_sigma_from_pseudobands.py &> script_update_sigma_from_pseudobands.py.out' if jmespath.search('parabands.enabled', self.wfn_options) else ''}

ln -sf {kgrid_link_folder} ./kgrid
ln -sf ../epsilon/epsmat.h5 ./
ln -sf ../epsilon/eps0mat.h5 ./
ln -sf ../{jmespath.search('gw.sigma.wfnlink', self.inputdict)}/rho ./RHO
ln -sf ../{jmespath.search('gw.sigma.wfnlink', self.inputdict)}/vxc.dat ./
ln -sf {wfnlink_str} ./WFN_inner.h5 
{scheduler.get_exec_prefix()}sigma.cplx.x &> sigma.inp.out
{unfold_string if jmespath.search('sym', self.wfn_options)==True else ""}
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        output = {
            './sigma/sigma.inp': self.sigma,
            './sigma/job_sigma.sh': self.job_sigma
        }
        
        # Add the script to from pseudobands if needed.
        if jmespath.search('parabands.enabled', self.wfn_options):
            output['./sigma/script_update_sigma_from_pseudobands.py'] = self.sigma_update_from_pseudobands

        return output
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './sigma/job_sigma.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './sigma',
        ]
#endregion