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
    def epsilon(self) -> str:
        # Qestruct.
        max_val_bands: int = QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('wfn.kgrid[0]', self.inputdict),
                jmespath.search('wfn.kgrid[1]', self.inputdict),
                jmespath.search('wfn.kgrid[2]', self.inputdict),
            ],
            qshift=[
                jmespath.search('wfnq.qshift[0]', self.inputdict),
                jmespath.search('wfnq.qshift[1]', self.inputdict),
                jmespath.search('wfnq.qshift[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfn.sym', self.inputdict),
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

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf {jmespath.search('gw.epsilon.wfnlink', self.inputdict)} ./WFN.h5 
ln -sf {jmespath.search('gw.epsilon.wfnqlink', self.inputdict)} ./WFNq.h5 
{scheduler.get_exec_prefix()}epsilon.cplx.x &> epsilon.inp.out 
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'epsilon.inp': self.epsilon,
            'job_epsilon.sh': self.job_epsilon,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_epsilon.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './epsilon.inp',
            './job_epsilon.sh',
            './WFN.h5',
            './WFNq.h5',
            './epsmat.h5',
            './eps0mat.h5',
            './epsilon.log',
            './chi_converge.dat',
            './epsilon.inp.out',
            './checkbz.log',
        ]
#endregion