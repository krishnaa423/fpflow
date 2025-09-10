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
import numpy as np 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeWfnStep(Step):
    @property
    def wfn(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )
        cond_bands: int = jmespath.search('wfn.cond_bands', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('wfn.kgrid[0]', self.inputdict),
                jmespath.search('wfn.kgrid[1]', self.inputdict),
                jmespath.search('wfn.kgrid[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfn.sym', self.inputdict),
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
        update_dict(wfndict, jmespath.search('wfn.args', self.inputdict))

        return QeGrammar().write(wfndict)
    
    @property
    def job_wfn(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wfn.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < wfn.in &> wfn.in.out 

cp ./tmp/struct.xml ./wfn.xml
'''
        return file_string
    
    @property
    def wfn_pw2bgw(self) -> str:
        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'WFN_coo'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('wfn.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('wfn.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('wfn.kgrid[2]', self.inputdict),
                'wfng_dk1': 0.0,
                'wfng_dk2': 0.0,
                'wfng_dk3': 0.0,
                'rhog_flag': '.true.',
                'rhog_file': "'RHO'",
                'vxcg_flag': '.true.',
                'vxcg_file': "'VXC'",
                'vscg_flag': '.true.',
                'vscg_file': "'VSC'",
                'vkbg_flag': '.true.',
                'vkbg_file': "'VKB'",
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('wfn.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)
    
    @property
    def job_wfn_pw2bgw(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wfn.job_pw2bgw_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw2bgw.x -pd .true. < wfn_pw2bgw.in &> wfn_pw2bgw.in.out
cp ./tmp/WFN_coo ./
cp ./tmp/RHO ./
cp ./tmp/VXC ./
cp ./tmp/VSC ./
cp ./tmp/VKB ./ 
'''
        return file_string
    
    @property
    def parabands(self) -> str:
        max_val_bands: int = QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )
        parabands_cond_bands: int = jmespath.search('wfn.parabands_cond_bands', self.inputdict)

        parabandsdict: dict = {
            'input_wfn_file': 'WFN_coo',
            'output_wfn_file': 'WFN_parabands.h5',
            'vsc_file': 'VSC',
            'vkb_file': 'VKB',
            'number_bands': parabands_cond_bands + max_val_bands,
            'wfn_io_mpiio_mode': 1,
        }

        # Update if needed. 
        update_dict(parabandsdict, jmespath.search('wfn.parabands_args', self.inputdict))

        return BgwGrammar().write(parabandsdict)
    
    @property
    def job_parabands(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wfn.job_parabands_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}parabands.cplx.x &> parabands.inp.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'wfn.in': self.wfn,
            'job_wfn.sh': self.job_wfn,
            'wfn_pw2bgw.in': self.wfn_pw2bgw,
            'job_wfn_pw2bgw.sh': self.job_wfn_pw2bgw,
            'parabands.inp': self.parabands,
            'job_parabands.sh': self.job_parabands,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_wfn.sh',
            './job_wfn_pw2bgw.sh',
            './job_parabands.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './wfn.in',
            './job_wfn.sh',
            './wfn_pw2bgw.in',
            './job_wfn_pw2bgw.sh',
            './parabands.inp',
            './job_parabands.sh',
            './tmp',
            './wfn.xml',
            './WFN_parabands.h5',
            './WFN_coo',
            './RHO',
            './VXC',
            './VSC',
            './VKB',
            './wfn.in.out',
            './wfn_pw2bgw.in.out',
            './parabands.inp.out',
        ]

#endregion