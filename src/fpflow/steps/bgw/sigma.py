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
    def sigma(self) -> str:
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
            is_reduced=jmespath.search('wfn.sym', self.inputdict),
        )

        sigmadict: dict = {
            'no_symmetries_q_grid': '',
            'number_bands': jmespath.search('gw.sigma.conv_cond_bands', self.inputdict) + max_val_bands,
            'band_index_min': max_val_bands - jmespath.search('gw.sigma.cond_bands', self.inputdict) + 1,
            'band_index_max': max_val_bands + jmespath.search('gw.sigma.val_bands', self.inputdict),
            'degeneracy_check_override': '',
            'screened_coulomb_cutoff': jmespath.search('gw.sigma.ecut', self.inputdict),
            'dont_use_vxcdat': '',
            'use_wfn_hdf5': '',
            'kpoints': kpts.sigma_kpts,
        }

        # Update if needed. 
        update_dict(sigmadict, jmespath.search('gw.epsilon.args', self.inputdict))

        return BgwGrammar().write(sigmadict)
    
    @property
    def job_sigma(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'gw.epsilon.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf {jmespath.search('gw.sigma.wfnlink', self.inputdict)} ./WFN_inner.h5 
{scheduler.get_exec_prefix()}sigma.cplx.x &> sigma.inp.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'sigma.inp': self.sigma,
            'job_sigma.sh': self.job_sigma
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_sigma.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './sigma.inp',
            './job_sigma.sh',
            './WFN_inner.h5',
            './eqp0.dat',
            './eqp1.dat',
            './sigma_hp.log',
            './ch_converge.dat',
            './sigma.inp.out',
            './dtmat',
            './vxc.dat',
            './x.dat',
        ]
#endregion