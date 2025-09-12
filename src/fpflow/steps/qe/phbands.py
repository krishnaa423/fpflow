#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.inputs.grammars.namelist import NamelistGrammar
import jmespath
from fpflow.io.update import update_dict
from fpflow.io.logging import get_logger
from fpflow.schedulers.scheduler import Scheduler
from importlib.util import find_spec
from fpflow.plots.phbands import PhbandsPlot
from fpflow.structure.kpath import Kpath
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QePhbandsStep(Step):
    @property
    def q2r_bands(self) -> str:
        q2rdict: dict = {
            'input': {
                'zasr': "'crystal'",
                'fildyn': "'struct.dyn'",
                'flfrc': "'struct.fc'",
            }
        }

        # Update if needed. 
        update_dict(q2rdict, jmespath.search('phbands.q2r_args', self.inputdict))

        return NamelistGrammar().write(q2rdict)
    
    @property
    def job_q2r_bands(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'phbands.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}q2r.x < q2r_bands.in &> q2r_bands.in.out
'''
        return file_string
    
    @property
    def matdyn_bands(self) -> str:
        matdyndict: dict = {
            'input': {
                'asr': "'crystal'",
                'flfrc': "'struct.fc'",
                'flfrq': "'struct.freq'",
                'flvec': "'struct.modes'",
                'q_in_band_form': '.true.',
                'q_in_cryst_coord': '.true.',
            }
        }

        # Update if needed. 
        update_dict(matdyndict, jmespath.search('phbands.matdyn_args', self.inputdict))

        matdyn_filestring: str =  NamelistGrammar().write(matdyndict)

        matdyn_filestring += Kpath.from_yamlfile().matdyn_str

        return matdyn_filestring
    
    @property
    def job_matdyn_bands(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'phbands.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}matdyn.x < matdyn_bands.in &> matdyn_bands.in.out 
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'q2r_bands.in': self.q2r_bands,
            'job_q2r_bands.sh': self.job_q2r_bands,
            'matdyn_bands.in': self.matdyn_bands,
            'job_matdyn_bands.sh': self.job_matdyn_bands,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_q2r_bands.sh',
            './job_matdyn_bands.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './q2r_bands.in',
            './q2r_bands.in.out',
            './job_q2r_bands.sh',
            './matdyn_bands.in',
            './matdyn_bands.in.out',
            './job_matdyn_bands.sh',
            './struct.dyn*',
            './struct.fc',
            './struct.freq',
            './struct.freq.gp',
            './struct.modes',
            'plot_phbands.h5'
        ]
    
    def plot(self, **kwargs):
        PhbandsPlot().save_plot(**kwargs)

#endregion