#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.inputs.grammars.namelist import NamelistGrammar
from fpflow.inputs.grammars.qe import QeGrammar
from fpflow.inputs.grammars.bgw import BgwGrammar
import jmespath
from fpflow.io.update import update_dict
from fpflow.io.logging import get_logger
from fpflow.schedulers.scheduler import Scheduler
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.structure.kpts import Kpts
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeWfnfiStep(Step):
    @property
    def wfnfi(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )
        cond_bands: int = jmespath.search('wfnfi.cond_bands', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid=[
                jmespath.search('wfnfi.kgrid[0]', self.inputdict),
                jmespath.search('wfnfi.kgrid[1]', self.inputdict),
                jmespath.search('wfnfi.kgrid[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfnfi.sym', self.inputdict),
        )

        wfndict: dict = {
            'control': {
                'outdir': './tmp',
                'prefix': 'struct',
                'pseudo_dir': './pseudos/qe',
                'calculation': 'bands',
                'tprnfor': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                'nbnd': int(cond_bands + max_val_bands),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict)
            },
            'electrons': {},
            'ions': {},
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'crystal',
                'nkpt': kpts.nkpt,
                'data': kpts.wfn_kpts,
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            wfndict['system']['noncolin'] = True
            wfndict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(wfndict, jmespath.search('wfnfi.args', self.inputdict))

        return QeGrammar().write(wfndict)

    @property
    def job_wfnfi(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wfnfi.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < wfnfi.in &> wfnfi.in.out
'''
        return file_string

    @property
    def wfnfi_pw2bgw(self) -> str:
        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'WFN_fii'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('wfnfi.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('wfnfi.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('wfnfi.kgrid[2]', self.inputdict),
                'wfng_dk1': 0.0,
                'wfng_dk2': 0.0,
                'wfng_dk3': 0.0,
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('wfnfi.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def job_wfnfi_pw2bgw(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wfnfi.job_pw2bgw_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw2bgw.x -pd .true. < wfnfi_pw2bgw.in &> wfnfi_pw2bgw.in.out
cp ./tmp/WFN_fii ./
cp ./tmp/struct.xml ./wfnfi.xml
wfn2hdf.x BIN WFN_fii WFN_fii.h5
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'wfnfi.in': self.wfnfi,
            'job_wfnfi.sh': self.job_wfnfi,
            'wfnfi_pw2bgw.in': self.wfnfi_pw2bgw,
            'job_wfnfi_pw2bgw.sh': self.job_wfnfi_pw2bgw,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_wfnfi.sh',
            './job_wfnfi_pw2bgw.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './wfnfi.in',
            './job_wfnfi.sh',
            './wfnfi_pw2bgw.in',
            './job_wfnfi_pw2bgw.sh',
            './tmp',
            './wfnfi.xml',
            './WFN_fii',
            './WFN_fii.h5',
            './wfnfi.in.out',
            './wfnfi_pw2bgw.in.out',
        ]

#endregion