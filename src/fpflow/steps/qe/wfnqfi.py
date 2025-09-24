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
class QeWfnqfiStep(Step):
    @property
    def wfnqfi(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )
        cond_bands: int = jmespath.search('wfnqfi.cond_bands', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid=[
                jmespath.search('wfnqfi.kgrid[0]', self.inputdict),
                jmespath.search('wfnqfi.kgrid[1]', self.inputdict),
                jmespath.search('wfnqfi.kgrid[2]', self.inputdict),
            ],
            qshift=[
                jmespath.search('wfnq.qshift[0]', self.inputdict),
                jmespath.search('wfnq.qshift[1]', self.inputdict),
                jmespath.search('wfnq.qshift[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfnqfi.sym', self.inputdict),
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
                'data': kpts.wfnq_kpts,
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            wfndict['system']['noncolin'] = True
            wfndict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(wfndict, jmespath.search('wfnqfi.args', self.inputdict))

        return QeGrammar().write(wfndict)

    @property
    def job_wfnqfi(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wfnqfi.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < wfnqfi.in &> wfnqfi.in.out 
'''
        return file_string

    @property
    def wfnqfi_pw2bgw(self) -> str:
        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'WFNq_fii'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('wfnqfi.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('wfnqfi.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('wfnqfi.kgrid[2]', self.inputdict),
                'wfng_dk1': jmespath.search('wfnq.qshift[0]', self.inputdict) * jmespath.search('wfnqfi.kgrid[0]', self.inputdict),
                'wfng_dk2': jmespath.search('wfnq.qshift[1]', self.inputdict) * jmespath.search('wfnqfi.kgrid[1]', self.inputdict),
                'wfng_dk3': jmespath.search('wfnq.qshift[2]', self.inputdict) * jmespath.search('wfnqfi.kgrid[2]', self.inputdict),
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('wfnqfi.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def job_wfnqfi_pw2bgw(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wfnqfi.job_pw2bgw_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw2bgw.x -pd .true. < wfnqfi_pw2bgw.in &> wfnqfi_pw2bgw.in.out
cp ./tmp/WFNq_fii ./
cp ./tmp/struct.xml ./wfnqfi.xml
wfn2hdf.x BIN WFNq_fii WFNq_fii.h5
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'wfnqfi.in': self.wfnqfi,
            'job_wfnqfi.sh': self.job_wfnqfi,
            'wfnqfi_pw2bgw.in': self.wfnqfi_pw2bgw,
            'job_wfnqfi_pw2bgw.sh': self.job_wfnqfi_pw2bgw,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_wfnqfi.sh',
            './job_wfnqfi_pw2bgw.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './wfnqfi.in',
            './job_wfnqfi.sh',
            './wfnqfi_pw2bgw.in',
            './job_wfnqfi_pw2bgw.sh',
            './tmp',
            './wfnqfi.xml',
            './WFNq_fii',
            './WFNq_fii.h5',
            './wfnqfi.in.out',
            './wfnqfi_pw2bgw.in.out',
        ]

#endregion