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
class BgwAbsorptionStep(Step):
    @property
    def absorption(self) -> str:
        qshift: list[int] = jmespath.search('bse.absorption.qshift[*]', self.inputdict)
        pol_dir: list[int] = jmespath.search('bse.absorption.pol_dir[*]', self.inputdict)

        absorptiondict: dict = {
            'exciton_Q_shift': f"2 {qshift[0]} {qshift[1]} {qshift[2]}",
            'use_symmetries_coarse_grid': '',
            'use_symmetries_fine_grid': '',
            'use_symmetries_shifted_grid': '',
            'number_val_bands_coarse': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_val_bands_fine': jmespath.search('bse.absorption.val_bands', self.inputdict) - 1,
            'number_cond_bands_coarse': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'number_cond_bands_fine': jmespath.search('bse.absorption.cond_bands', self.inputdict),
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
    
    @property
    def job_absorption(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bse.absorption.job_info')
        link_dir_prefix: str = jmespath.search('bse.absorption.link_dir_prefix', self.inputdict)

        if link_dir_prefix is not None:
            file_string = f'''#!/bin/bash
{scheduler.get_script_header()}


ln -sf {os.path.join(link_dir_prefix, 'epsmat.h5')} ./epsmat.h5 
ln -sf {os.path.join(link_dir_prefix, 'eps0mat.h5')} ./eps0mat.h5 
ln -sf {os.path.join(link_dir_prefix, 'eqp1.dat')} ./eqp_co.dat 
ln -sf {os.path.join(link_dir_prefix, jmespath.search('bse.absorption.wfnco_link', self.inputdict))} ./WFN_co.h5 
ln -sf {os.path.join(link_dir_prefix, jmespath.search('bse.absorption.wfnqco_link', self.inputdict))} ./WFNq_co.h5 
ln -sf {os.path.join(link_dir_prefix, jmespath.search('bse.absorption.wfnfi_link', self.inputdict))} ./WFN_fi.h5 
ln -sf {os.path.join(link_dir_prefix, jmespath.search('bse.absorption.wfnqfi_link', self.inputdict))} ./WFNq_fi.h5 
{scheduler.get_exec_prefix()}absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
'''
        else:
            file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf {jmespath.search('bse.absorption.wfnco_link', self.inputdict)} ./WFN_co.h5 
ln -sf {jmespath.search('bse.absorption.wfnqco_link', self.inputdict)} ./WFNq_co.h5 
ln -sf {jmespath.search('bse.absorption.wfnfi_link', self.inputdict)} ./WFN_fi.h5 
ln -sf {jmespath.search('bse.absorption.wfnqfi_link', self.inputdict)} ./WFNq_fi.h5 
ln -sf eqp1.dat eqp_co.dat 
{scheduler.get_exec_prefix()}absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
'''
        return file_string
    
    @property
    def file_contents(self) -> dict:
        return {
            'absorption.inp': self.absorption,
            'job_absorption.sh': self.job_absorption,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_absorption.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './absorption.inp',
            './job_absorption.sh',
            './WFN_co.h5',
            './WFNq_co.h5',
            './WFN_fi.h5',
            './WFNq_fi.h5',
            './eigenvalues.dat',
            './eigenvalues_noeh.dat',
            './absorption_eh.dat',
            './absorption_noeh.dat',
            './dvmat_norm.dat',
            './dcmat_norm.dat',
            './eqp_co.dat',
            './eqp.dat',
            './eqp_q.dat',
            './bandstructure_absorption.dat',
            './eigenvectors.h5',
            './hbse*.h5',
            './x.dat',
            './epsdiag.dat',
            './dtmat',
            './vmtxel',
            './absorption.inp.out',
        ]
#endregion