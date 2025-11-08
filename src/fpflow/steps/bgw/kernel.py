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
            'exciton_Q_shift': f"2 {qshift[0]:15.10f} {qshift[1]:15.10f} {qshift[2]:15.10f}",
            'use_symmetries_coarse_grid': '',
            'number_val_bands': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_cond_bands': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'use_wfn_hdf5': '',
            'dont_check_norms': '',
            'energy_loss': '',
        }

        # Update if needed. 
        update_dict(kerneldict, jmespath.search('bse.kernel.args', self.inputdict))

        return BgwGrammar().write(kerneldict)

    def get_link(self, path: str) -> str:
        wfnlink: str = jmespath.search(path, self.inputdict)
        wfn_names: list[str] = jmespath.search(f'nscf.list[*].name', self.inputdict)
        wfn_index: int = wfn_names.index(wfnlink)
        self.wfn_options: dict = jmespath.search(f'nscf.list[{wfn_index}]', self.inputdict)
        wfnlink_str: str = os.path.join(
            jmespath.search('bse.link_dir_prefix', self.inputdict),
            jmespath.search(path, self.inputdict),
            'wfn_parabands.h5' if jmespath.search('parabands.enabled', self.wfn_options) else 'wfn.h5' 
        )

        return wfnlink_str

    @property
    def job_kernel(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bse.kernel.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf {jmespath.search('bse.link_dir_prefix', self.inputdict)}/epsilon/epsmat.h5 ./
ln -sf {jmespath.search('bse.link_dir_prefix', self.inputdict)}/epsilon/eps0mat.h5 ./
ln -sf {self.get_link('bse.absorption.wfnco_link')} ./WFN_co.h5 
ln -sf {self.get_link('bse.absorption.wfnqco_link')} ./WFNq_co.h5 
{scheduler.get_exec_prefix()}kernel.cplx.x &> kernel.inp.out
    '''

        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './kernel/kernel.inp': self.kernel,
            './kernel/job_kernel.sh': self.job_kernel,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './kernel/job_kernel.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './kernel',
        ]
#endregion