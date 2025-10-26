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
class BgwEpsilonStep(Step):
    @property
    def eps_update_from_pseudobands(self) -> str:
        return f'''#!/usr/bin/env python

import h5py
import glom 
from fpflow.inputs.inputyaml import InputYaml
from fpflow.steps.bgw.epsilon import BgwEpsilonStep
from fpflow.structure.qe.qe_struct import QeStruct
import jmespath
import os

inputdict: dict = InputYaml.from_yaml_file('../input.yaml').inputdict
num_pseudobands: int = 0
wfnlink: str = jmespath.search('gw.epsilon.wfnlink', inputdict)
parabands_file: str = os.path.join('..', wfnlink, 'wfn_parabands.h5')
with h5py.File(parabands_file, 'r') as f:
    num_pseudobands = f['/parabands/pseudobands/nb_total'][()]
# Qestruct.
max_val_bands: int = int(QeStruct.from_inputdict(inputdict).max_val(
    xc=jmespath.search('scf.xc', inputdict),
    is_soc=jmespath.search('scf.is_spinorbit', inputdict),
))

num_gw_cond_bands: int = num_pseudobands - max_val_bands - 1
glom.assign(inputdict, 'gw.epsilon.cond_bands', num_gw_cond_bands)
eps_step = BgwEpsilonStep(inputdict=inputdict)
os.chdir('../')
eps_step.create()
os.chdir('./epsilon')
'''

    @property
    def epsilon(self) -> str:
        # Qestruct.
        max_val_bands: int = int(QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))

        # wfnlink. 
        wfn_link: str = jmespath.search('gw.epsilon.wfnlink', self.inputdict)
        wfn_names: list[str] = jmespath.search(f'nscf.list[*].name', self.inputdict)
        wfn_index: int = wfn_names.index(wfn_link)
        self.wfn_options: dict = jmespath.search(f'nscf.list[{wfn_index}]', self.inputdict)

        # wfnqlink.
        wfnq_link: str = jmespath.search('gw.epsilon.wfnqlink', self.inputdict)
        wfnq_names: list[str] = jmespath.search(f'nscf.list[*].name', self.inputdict)
        wfnq_index: int = wfnq_names.index(wfnq_link)
        self.wfnq_options: dict = jmespath.search(f'nscf.list[{wfnq_index}]', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('kgrid[0]', self.wfnq_options),
                jmespath.search('kgrid[1]', self.wfnq_options),
                jmespath.search('kgrid[2]', self.wfnq_options),
            ],
            qshift=[
                jmespath.search('qshift[0]', self.wfnq_options),
                jmespath.search('qshift[1]', self.wfnq_options),
                jmespath.search('qshift[2]', self.wfnq_options),
            ],
            is_reduced=jmespath.search('sym', self.wfnq_options),
        )

        epsilondict: dict = {
            'number_bands': jmespath.search('gw.epsilon.cond_bands', self.inputdict) + max_val_bands,
            'degeneracy_check_override': '',
            'epsilon_cutoff': jmespath.search('gw.epsilon.ecut', self.inputdict),
            'use_wfn_hdf5': '',
            'qpoints': kpts.epsilon_kpts 
        }

        # Update if needed. 
        update_dict(epsilondict, jmespath.search('gw.epsilon.args', self.inputdict))

        return BgwGrammar().write(epsilondict)
    
    @property
    def job_epsilon(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'gw.epsilon.job_info')

        wfnlink_str: str = os.path.join(
            '..',
            jmespath.search('gw.epsilon.wfnlink', self.inputdict),
            'wfn_parabands.h5' if jmespath.search('parabands.enabled', self.wfn_options) else 'wfn.h5' 
        )

        wfnqlink_str: str = os.path.join(
            '..',
            jmespath.search('gw.epsilon.wfnqlink', self.inputdict),
            'wfn_parabands.h5' if jmespath.search('parabands.enabled', self.wfnq_options) else 'wfn.h5' 
        )

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{'python ./script_update_eps_from_pseudobands.py &> script_update_eps_from_pseudobands.py.out' if jmespath.search('parabands.enabled', self.wfn_options) else ''}

ln -sf {wfnlink_str} ./WFN.h5 
ln -sf {wfnqlink_str} ./WFNq.h5 
{scheduler.get_exec_prefix()}epsilon.cplx.x &> epsilon.inp.out 
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        output = {
            './epsilon/epsilon.inp': self.epsilon,
            './epsilon/job_epsilon.sh': self.job_epsilon,
        }

        # Add the script to from pseudobands if needed.
        if jmespath.search('parabands.enabled', self.wfn_options):
            output['./epsilon/script_update_eps_from_pseudobands.py'] = self.eps_update_from_pseudobands

        return output
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './epsilon/job_epsilon.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './epsilon',
        ]
#endregion