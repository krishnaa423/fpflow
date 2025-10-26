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
from fpflow.plots.kpdos import KpdosPlot 

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeKpdosStep(Step):
    @property
    def kpdos(self):
        kpdosdict: dict = {
            'projwfc': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'kresolveddos': '.true.',
                'filpdos': "'kpdos.dat'",
            }
        }

        # Update if needed. 
        update_dict(kpdosdict, jmespath.search('kpdos.args', self.inputdict))

        file_string: str =  NamelistGrammar().write(kpdosdict)

        return file_string

    @property
    def job_kpdos(self):
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'kpdos.job_info')
        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf ../dftelbands/tmp ./tmp
{scheduler.get_exec_prefix()}projwfc.x -pd .true. {scheduler.get_exec_infix()} < kpdos.in &> kpdos.in.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './kpdos/kpdos.in': self.kpdos,
            './kpdos/job_kpdos.sh': self.job_kpdos,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './kpdos/job_kpdos.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './kpdos',
        ]
    
    def plot(self, **kwargs):
            KpdosPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion