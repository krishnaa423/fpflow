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
class BgwKernelStep(Step):
    @property
    def kernel(self) -> str:
        qshift: list[int] = jmespath.search('bse.kernel.qshift[*]', self.inputdict)

        kerneldict: dict = {
            'exciton_Q_shift': f"2 {qshift[0]} {qshift[1]} {qshift[2]}",
            'use_symmetries_coarse_grid': '',
            'number_val_bands': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_cond_bands': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'use_wfn_hdf5': '',
            'dont_check_norms': '',
        }

        # Update if needed. 
        update_dict(kerneldict, jmespath.search('bse.kernel.args', self.inputdict))

        return BgwGrammar().write(kerneldict)

    @property
    def job_kernel(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bse.kernel.job_info')
        link_dir_prefix: str = jmespath.search('bse.kernel.link_dir_prefix', self.inputdict)

        if link_dir_prefix is not None:
            file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf {os.path.join(link_dir_prefix, 'epsmat.h5')} ./epsmat.h5
ln -sf {os.path.join(link_dir_prefix, 'eps0mat.h5')} ./eps0mat.h5
ln -sf {os.path.join(link_dir_prefix, jmespath.search('bse.absorption.wfnco_link', self.inputdict))} ./WFN_co.h5 
ln -sf {os.path.join(link_dir_prefix, jmespath.search('bse.absorption.wfnqco_link', self.inputdict))} ./WFNq_co.h5 
{scheduler.get_exec_prefix()}kernel.cplx.x &> kernel.inp.out
    '''
        else:
            file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf {jmespath.search('bse.absorption.wfnco_link', self.inputdict)} ./WFN_co.h5 
ln -sf {jmespath.search('bse.absorption.wfnqco_link', self.inputdict)} ./WFNq_co.h5 
{scheduler.get_exec_prefix()}kernel.cplx.x &> kernel.inp.out
    '''

        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'kernel.inp': self.kernel,
            'job_kernel.sh': self.job_kernel,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_kernel.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './kernel.inp',
            './job_kernel.sh',
            './WFN_co.h5',
            './WFNq_co.h5',
            './bsemat.h5',
            './kernel.inp.out',
        ]
#endregion