#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.schedulers.scheduler import Scheduler
from fpflow.steps.step import Step 
from fpflow.io.update import update_dict
import jmespath
from fpflow.inputs.grammars.namelist import NamelistGrammar
from fpflow.structure.kpath import Kpath
from fpflow.inputs.grammars.qe import QeGrammar
from fpflow.structure.kpts import Kpts
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.plots.pdos import PdosPlot

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QePdosStep(Step):
    @property
    def pdos(self):
        pdosdict: dict = {
            'projwfc': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'filpdos': "'pdos.dat'",
            }
        }

        # Update if needed. 
        update_dict(pdosdict, jmespath.search('pdos.args', self.inputdict))

        file_string: str =  NamelistGrammar().write(pdosdict)

        return file_string

    @property
    def job_pdos(self):
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'pdos.job_info')
        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf ../dos/tmp ./tmp
{scheduler.get_exec_prefix()}projwfc.x -pd .true. {scheduler.get_exec_infix()} < pdos.in &> pdos.in.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './pdos/pdos.in': self.pdos,
            './pdos/job_pdos.sh': self.job_pdos,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './pdos/job_pdos.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
             './pdos',
        ]
    
    def plot(self, **kwargs):
            PdosPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion