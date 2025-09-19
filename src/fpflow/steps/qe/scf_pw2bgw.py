#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.schedulers.scheduler import Scheduler
import jmespath
from fpflow.io.update import update_dict
from fpflow.inputs.grammars.namelist import NamelistGrammar

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeScfPw2bgwStep(Step):
    @property
    def pw2bgw(self) -> str:
        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'WFN_scf'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('scf.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('scf.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('scf.kgrid[2]', self.inputdict),
                'wfng_dk1': 0.0,
                'wfng_dk2': 0.0,
                'wfng_dk3': 0.0,
                'rhog_flag': '.true.',
                'rhog_file': "'RHO_scf'",
                'vxcg_flag': '.true.',
                'vxcg_file': "'VXCG_scf'",
                'vscg_flag': '.true.',
                'vscg_file': "'VSCG_scf'",
                'vkbg_flag': '.true.',
                'vkbg_file': "'VKBG_scf'",
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('scf_pw2bgw.args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def job_pw2bgw(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'scf_pw2bgw.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw2bgw.x -pd .true. < scf_pw2bgw.in &> scf_pw2bgw.in.out
cp ./tmp/WFN_scf ./WFN_scf
cp ./tmp/RHO_scf ./
cp ./tmp/VXC_scf ./
cp ./tmp/VSC_scf ./
cp ./tmp/VKB_scf ./
wfn2hdf.x BIN WFN_scf WFN_scf.h5 
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'scf_pw2bgw.in': self.pw2bgw,
            'job_scf_pw2bgw.sh': self.job_pw2bgw,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_scf_pw2bgw.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './scf_pw2bgw.in',
            './job_scf_pw2bgw.sh',
            './WFN_scf',
            './WFN_scf.h5',
            './RHO_scf',
            './VXC_scf',
            './VSC_scf',
            './VKB_scf',
        ]

#endregion