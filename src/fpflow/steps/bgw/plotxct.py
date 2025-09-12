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
class BgwPlotxctStep(Step):
    @property
    def plotxct(self) -> str:
        qshift: list[int] = jmespath.search('bse.absorption.qshift[*]', self.inputdict)

        plotxctdict: dict = {
            'hole_position': ' '.join(list(map(str, jmespath.search('bse.plotxct.hole_position', self.inputdict)))),
            'supercell_size': ' '.join(list(map(str, jmespath.search('bse.plotxct.supercell_size', self.inputdict)))),
            'use_symmetries_fine_grid': '',
            'use_symmetries_shifted_grid': '',
            'plot_spin': 1,
            'plot_state': jmespath.search('bse.plotxct.xct_state', self.inputdict),
            'use_wfn_hdf5': '',
        }

        # Update if needed. 
        update_dict(plotxctdict, jmespath.search('bse.plotxct.args', self.inputdict))

        return BgwGrammar().write(plotxctdict)
    
    @property
    def job_plotxct(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bse.plotxct.job_info')
        link_dir_prefix: str = jmespath.search('bse.absorption.link_dir_prefix', self.inputdict)

        if link_dir_prefix is not None:
            file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf {os.path.join(link_dir_prefix, jmespath.search('bse.absorption.wfnfi_link', self.inputdict))} ./WFN_fi.h5 
ln -sf {os.path.join(link_dir_prefix, jmespath.search('bse.absorption.wfnqfi_link', self.inputdict))} ./WFNq_fi.h5 
{scheduler.get_exec_prefix()}plotxct.cplx.x &> plotxct.inp.out 
volume.py {os.path.join(link_dir_prefix, './scf.in')} espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
rm -rf *.a3Dr
'''
        else:
            file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf {jmespath.search('bse.absorption.wfnfi_link', self.inputdict)} ./WFN_fi.h5 
ln -sf {jmespath.search('bse.absorption.wfnqfi_link', self.inputdict)} ./WFNq_fi.h5 
{scheduler.get_exec_prefix()}plotxct.cplx.x &> plotxct.inp.out 
volume.py ./scf.in espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
rm -rf *.a3Dr
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'plotxct.inp': self.plotxct,
            'job_plotxct.sh': self.job_plotxct,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_plotxct.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './plotxct.inp',
            './job_plotxct.sh',
            './*.a3Dr',
            './plotxcf.xsf',
            './plotxct_elec*.xsf',
            './plotxct_hole*.xsf',
            './plotxct.inp.out',
            './WFN_fi.h5',
            './WFNq_fi.h5',
        ]
#endregion