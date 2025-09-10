#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.schedulers.scheduler import Scheduler
from fpflow.plots.xctph import XctphPlot
from fpflow.structure.qe.qe_struct import QeStruct
import jmespath
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class XctphStep(Step):
    @property
    def xctph(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        nocc: int = int(qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))
        nc: int = jmespath.search('bse.absorption.cond_bands', self.inputdict)
        nv: int = jmespath.search('bse.absorption.val_bands', self.inputdict)
        nxct: int = jmespath.search('bse.absorption.nxct', self.inputdict)
        npool: int = 1 if jmespath.search('epw.job_info.nk', self.inputdict) is None else jmespath.search('epw.job_info.nk', self.inputdict)


        file_string: str = f'''
#region modules
from xctph.xctph import Xctph
#endregion

#region variables
#endregion

#region functions
def main():
    xctph: Xctph = Xctph(nocc={nocc}, nc={nc}, nv={nv}, nxct={nxct}, npool={npool})
    xctph.calc()
    xctph.write()
#endregion

#region classes
#endregion

#region main
main()
#endregion
'''
        
        return file_string

    @property
    def job_xctph(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'xctph.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

rm -rf xctph.out
touch xctph.out
exec &> xctph.out

{scheduler.get_exec_prefix()}python3 script_xctph.py &> script_xctph.py.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'script_xctph.py': self.xctph,
            'job_xctph.sh': self.job_xctph, 
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_xctph.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './script_xctph.py',
            './job_xctph.sh',
            './xct.h5',
            './elph*.h5',
            './xctph*.h5',
            './xctph.out',
            './script_xctph.py.out',
        ] 
    
    def plot(self, **kwargs):
        XctphPlot().save_plot(**kwargs)

#endregion