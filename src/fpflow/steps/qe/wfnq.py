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
class QeWfnqStep(Step):
    @property
    def wfnq(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )
        cond_bands: int = jmespath.search('wfnq.cond_bands', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid=[
                jmespath.search('wfnq.kgrid[0]', self.inputdict),
                jmespath.search('wfnq.kgrid[1]', self.inputdict),
                jmespath.search('wfnq.kgrid[2]', self.inputdict),
            ],
            qshift=[
                jmespath.search('wfnq.qshift[0]', self.inputdict),
                jmespath.search('wfnq.qshift[1]', self.inputdict),
                jmespath.search('wfnq.qshift[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfnq.sym', self.inputdict),
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
        update_dict(wfndict, jmespath.search('wfn.args', self.inputdict))

        return QeGrammar().write(wfndict)

    @property
    def job_wfnq(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wfnq.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < wfnq.in &> wfnq.in.out 
'''
        return file_string

    @property
    def wfnq_pw2bgw(self) -> str:
        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'WFNq_coo'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('wfnq.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('wfnq.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('wfnq.kgrid[2]', self.inputdict),
                'wfng_dk1': jmespath.search('wfnq.qshift[0]', self.inputdict),
                'wfng_dk2': jmespath.search('wfnq.qshift[1]', self.inputdict),
                'wfng_dk3': jmespath.search('wfnq.qshift[2]', self.inputdict),
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('wfnq.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def job_wfnq_pw2bgw(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wfnq.job_pw2bgw_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw2bgw.x -pd .true. < wfnq_pw2bgw.in &> wfnq_pw2bgw.in.out
cp ./tmp/WFNq_coo ./
cp ./tmp/struct.xml ./wfnq.xml
wfn2hdf.x BIN WFNq_coo WFNq_coo.h5  
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'wfnq.in': self.wfnq,
            'job_wfnq.sh': self.job_wfnq,
            'wfnq_pw2bgw.in': self.wfnq_pw2bgw,
            'job_wfnq_pw2bgw.sh': self.job_wfnq_pw2bgw,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_wfnq.sh',
            './job_wfnq_pw2bgw.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './wfnq.in',
            './job_wfnq.sh',
            './wfnq_pw2bgw.in',
            './job_wfnq_pw2bgw.sh',
            './tmp',
            './wfnq.xml',
            './WFNq_coo',
            './WFNq_coo.h5',
            './wfnq.in.out',
            './wfnq_pw2bgw.in.out',
        ]
    
#endregion