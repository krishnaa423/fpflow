#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.inputs.grammars.qe import QeGrammar
from fpflow.inputs.grammars.namelist import NamelistGrammar
from fpflow.structure.qe.qe_struct import QeStruct
import jmespath
from fpflow.io.update import update_dict
from fpflow.io.logging import get_logger
from fpflow.schedulers.scheduler import Scheduler
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeScfStep(Step):
    @property
    def scf(self):
        qestruct = QeStruct.from_inputdict(self.inputdict)

        qedict: dict = {
            'control': {
                'outdir': './tmp',
                'prefix': 'struct',
                'pseudo_dir': './pseudos',
                'calculation': 'scf',
                'tprnfor': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict)
            },
            'electrons': {},
            'ions': {},
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'automatic',
                'data': [
                    jmespath.search('scf.kgrid[0]', self.inputdict),
                    jmespath.search('scf.kgrid[1]', self.inputdict),
                    jmespath.search('scf.kgrid[2]', self.inputdict),
                    0,
                    0,
                    0,
                ],
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            qedict['system']['noncolin'] = True
            qedict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(qedict, jmespath.search('scf.args', self.inputdict))

        return QeGrammar().write(qedict)

    @property
    def job_scf(self):
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'scf.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf ../pseudos/qe ./pseudos
{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < scf.in &> scf.in.out

cp ./tmp/struct.save/data-file-schema.xml ./scf.xml
'''
        return file_string

    @property
    def pw2bgw(self) -> str:
        # Qestruct.
        max_val_bands: int = QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )

        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'wfn'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('scf.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('scf.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('scf.kgrid[2]', self.inputdict),
                'wfng_dk1': 0.0,
                'wfng_dk2': 0.0,
                'wfng_dk3': 0.0,
                'rhog_flag': '.true.',
                'rhog_file': "'rho'",
                'vxc_flag': '.true.',
                'vxc_file': "'vxc.dat'",
                'vxc_diag_nmin': 1,
                'vxc_diag_nmax': max_val_bands,
                'vxc_offdiag_nmin': 0,
                'vxc_offdiag_nmax': 0,
                'vscg_flag': '.true.',
                'vscg_file': "'vsc'",
                'vkbg_flag': '.true.',
                'vkbg_file': "'vkb'",
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('scf.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def job_pw2bgw(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'scf.job_pw2bgw_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw2bgw.x -pd .true. < pw2bgw.in &> pw2bgw.in.out
cp ./tmp/wfn ./
cp ./tmp/rho ./
cp ./tmp/vxc.dat ./
cp ./tmp/vsc ./
cp ./tmp/vkb ./
wfn2hdf.x BIN wfn wfn.h5
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        contents_dict: dict = {
            './scf/scf.in': self.scf,
            './scf/job_scf.sh': self.job_scf,
        }

        # Add pw2bgw files if needed.
        if jmespath.search('scf.run_pw2bgw', self.inputdict):
            contents_dict.update({
                './scf/pw2bgw.in': self.pw2bgw,
                './scf/job_pw2bgw.sh': self.job_pw2bgw,
            })

        return contents_dict
    
    @property
    def job_scripts(self) -> List[str]:
        scrips_list: List[str] = [
            './scf/job_scf.sh'
        ]

        # Add pw2bgw script if needed.
        if jmespath.search('scf.run_pw2bgw', self.inputdict):
            scrips_list.append('./scf/job_pw2bgw.sh')

        return scrips_list

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './scf',
        ]

#endregion