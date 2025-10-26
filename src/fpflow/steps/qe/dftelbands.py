#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.inputs.grammars.qe import QeGrammar
from fpflow.structure.qe.qe_struct import QeStruct
import jmespath
from fpflow.io.update import update_dict
from fpflow.io.logging import get_logger
from fpflow.schedulers.scheduler import Scheduler
from fpflow.structure.kpath import Kpath
from fpflow.inputs.grammars.namelist import NamelistGrammar
from fpflow.plots.dftelbands import DftelbandsPlot

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeDftelbandsStep(Step):
    @property
    def dftelbands(self):
        #qestruct. 
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )
        cond_bands: int = jmespath.search('dftelbands.cond_bands', self.inputdict)

        # Kpath.
        kpath_data: list = Kpath.from_yamlfile().dftelbands_list

        qedict: dict = {
            'control': {
                'outdir': './tmp',
                'prefix': 'struct',
                'pseudo_dir': './pseudos',
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
                'type': 'crystal_b',
                'nkpt': len(kpath_data),
                'data': kpath_data,
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            qedict['system']['noncolin'] = True
            qedict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(qedict, jmespath.search('dftelbands.args', self.inputdict))

        return QeGrammar().write(qedict)

    @property
    def dftelbands_pw2bgw(self):
        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'wfn'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': 0,
                'wfng_nk2': 0,
                'wfng_nk3': 0,
                'wfng_dk1': 0.0,
                'wfng_dk2': 0.0,
                'wfng_dk3': 0.0,
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('dftelbands.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def job_dftelbands(self):
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'dftelbands.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

rm -rf ./tmp
cp -r ../scf/tmp ./tmp
{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < dftelbands.in &> dftelbands.in.out

cp ./tmp/struct.xml ./dftelbands.xml
'''
        return file_string

    @property
    def job_dftelbands_pw2bgw(self):
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'dftelbands.job_pw2bgw_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_infix()}pw2bgw.x -pd .true. < pw2bgw.in &> pw2bgw.in.out 
cp ./tmp/wfn ./
wfn2hdf.x BIN wfn wfn.h5
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './dftelbands/dftelbands.in': self.dftelbands,
            './dftelbands/pw2bgw.in': self.dftelbands_pw2bgw,
            './dftelbands/job_dftelbands.sh': self.job_dftelbands,
            './dftelbands/job_pw2bgw.sh': self.job_dftelbands_pw2bgw,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './dftelbands/job_dftelbands.sh',
            './dftelbands/job_pw2bgw.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './dftelbands',
        ]
    
    def plot(self, **kwargs):
        DftelbandsPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion