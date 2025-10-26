#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.schedulers.scheduler import Scheduler
from fpflow.inputs.grammars.bgw import BgwGrammar
import jmespath
from fpflow.io.update import update_dict
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.plots.gwelbands import GwelbandsPlot

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class BgwGwelbandsStep(Step):
    @property
    def inteqp(self):
        cond_bands: int = jmespath.search('gw.gwelbands.cond_bands', self.inputdict)
        val_bands: int = jmespath.search('gw.gwelbands.val_bands', self.inputdict)

        epsilondict: dict = {
            'number_val_bands_coarse': val_bands,
            'number_cond_bands_coarse': cond_bands,
            'number_val_bands_fine': val_bands,
            'number_cond_bands_fine': cond_bands,
            'degeneracy_check_override': '',
            'use_symmetries_coarse_grid': '',
            'no_symmetries_fine_grid': '',
        }

        # Update if needed. 
        update_dict(epsilondict, jmespath.search('gw.gwelbands.args', self.inputdict))

        return BgwGrammar().write(epsilondict)

    @property
    def job_inteqp(self):
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'gw.gwelbands.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf ../{jmespath.search('gw.gwelbands.wfnco_link', self.inputdict)}/wfn ./WFN_co 
ln -sf ../{jmespath.search('gw.gwelbands.wfnfi_link', self.inputdict)}/wfn ./WFN_fi 
ln -sf ../sigma/eqp1.dat ./eqp_co.dat 
{scheduler.get_exec_prefix()}inteqp.cplx.x &> inteqp.inp.out 
'''
        return file_string
    
    @property
    def file_contents(self) -> dict:
        return {
            './gwelbands/inteqp.inp': self.inteqp,
            './gwelbands/job_inteqp.sh': self.job_inteqp,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './gwelbands/job_inteqp.sh',   
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './gwelbands',
        ]
    
    def plot(self, **kwargs):
            GwelbandsPlot(inputdict=self.inputdict).save_figures(**kwargs)
    
#endregion