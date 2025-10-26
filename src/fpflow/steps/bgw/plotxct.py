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

        # Add spinor if needed. 
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            plotxctdict['spinor'] = ''
            plotxctdict['electron_spin'] = jmespath.search('bse.plotxct.spinor.electron_spin', self.inputdict)
            plotxctdict['hole_spin'] = jmespath.search('bse.plotxct.spinor.hole_spin', self.inputdict)

        # Update if needed. 
        update_dict(plotxctdict, jmespath.search('bse.plotxct.args', self.inputdict))

        return BgwGrammar().write(plotxctdict)
    
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
    def job_plotxct(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bse.plotxct.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf {self.get_link('bse.absorption.wfnfi_link')} ./WFN_fi.h5 
ln -sf {self.get_link('bse.absorption.wfnqfi_link')} ./WFNq_fi.h5 
ln -sf {jmespath.search('bse.link_dir_prefix', self.inputdict)}/absorption/eigenvectors.h5 ./
{scheduler.get_exec_prefix()}plotxct.cplx.x &> plotxct.inp.out 
volume.py {jmespath.search('bse.link_dir_prefix', self.inputdict)}/scf/scf.in espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
rm -rf *.a3Dr
'''
        
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './plotxct/plotxct.inp': self.plotxct,
            './plotxct/job_plotxct.sh': self.job_plotxct,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './plotxct/job_plotxct.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './plotxct',
        ]
#endregion