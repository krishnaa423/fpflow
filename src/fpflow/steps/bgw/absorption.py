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
from fpflow.plots.absorption import BseAbsorptionPlot
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class BgwAbsorptionStep(Step):
    @property
    def absorption(self) -> str:
        qshift: list[int] = jmespath.search('bse.absorption.qshift[*]', self.inputdict)
        pol_dir: list[int] = jmespath.search('bse.absorption.pol_dir[*]', self.inputdict)

        absorptiondict: dict = {
            'exciton_Q_shift': f"2 {qshift[0]:15.10f} {qshift[1]:15.10f} {qshift[2]:15.10f}",
            'use_symmetries_coarse_grid': '',
            'use_symmetries_fine_grid': '',
            'use_symmetries_shifted_grid': '',
            'number_val_bands_coarse': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_val_bands_fine': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_cond_bands_coarse': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'number_cond_bands_fine': jmespath.search('bse.absorption.cond_bands', self.inputdict) - 1,
            'degeneracy_check_override': '',
            'diagonalization': '',
            # 'use_elpa': '',  
            'use_momentum': '',  
            # 'use_velocity': '',  
            'polarization': ' '.join(list(map(str, pol_dir))),
            'eqp_co_corrections': '',
            'dump_bse_hamiltonian': '',
            'use_wfn_hdf5': '',
            'energy_resolution': 0.1,
            'write_eigenvectors': jmespath.search('bse.absorption.nxct', self.inputdict),
            'dont_check_norms': '',
        }

        # Update if needed. 
        update_dict(absorptiondict, jmespath.search('bse.absorption.args', self.inputdict))

        return BgwGrammar().write(absorptiondict)
    
    def get_link(self, path: str) -> str:
        wfnlink: str = jmespath.search(path, self.inputdict)
        wfn_names: list[str] = jmespath.search(f'nscf.list[*].name', self.inputdict)
        wfn_index: int = wfn_names.index(wfnlink)
        self.wfn_options: dict = jmespath.search(f'nscf.list[{wfn_index}]', self.inputdict)
        wfnlink_str: str = os.path.join(
            '..',
            jmespath.search(path, self.inputdict),
            'wfn_parabands.h5' if jmespath.search('parabands.enabled', self.wfn_options) else 'wfn.h5' 
        )

        return wfnlink_str

    @property
    def job_absorption(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bse.absorption.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}


ln -sf ../epsilon/epsmat.h5 ./
ln -sf ../epsilon/eps0mat.h5 ./
ln -sf ../sigma/eqp1.dat ./eqp_co.dat
ln -sf ../kernel/bsemat.h5 ./
ln -sf {self.get_link('bse.absorption.wfnco_link')} ./WFN_co.h5 
ln -sf {self.get_link('bse.absorption.wfnqco_link')} ./WFNq_co.h5 
ln -sf {self.get_link('bse.absorption.wfnfi_link')} ./WFN_fi.h5 
ln -sf {self.get_link('bse.absorption.wfnqfi_link')} ./WFNq_fi.h5 
{scheduler.get_exec_prefix()}absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
'''
        return file_string
    
    @property
    def file_contents(self) -> dict:
        return {
            './absorption/absorption.inp': self.absorption,
            './absorption/job_absorption.sh': self.job_absorption,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './absorption/job_absorption.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './absorption',
        ]
    
    def plot(self, **kwargs):
        BseAbsorptionPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion