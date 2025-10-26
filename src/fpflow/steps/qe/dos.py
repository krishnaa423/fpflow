#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.schedulers.scheduler import Scheduler
from fpflow.steps.step import Step 
from fpflow.io.update import update_dict
import jmespath
from fpflow.inputs.grammars.namelist import NamelistGrammar
from fpflow.structure.kpath import Kpath
from fpflow.inputs.grammars.qe import QeGrammar
from fpflow.structure.kpts import Kpts
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.plots.dos import DosPlot

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeDosStep(Step):
    @property
    def wfndos(self):
         # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = int(qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))
        cond_bands: int = jmespath.search('dos.cond_bands', self.inputdict)

        wfndosdict: dict = {
            'control': {
                'outdir': './tmp',
                'prefix': 'struct',
                'pseudo_dir': './pseudos',
                'calculation': 'bands',
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
                'type': 'automatic',
                'data': [
                    jmespath.search('dos.kgrid[0]', self.inputdict),
                    jmespath.search('dos.kgrid[1]', self.inputdict),
                    jmespath.search('dos.kgrid[2]', self.inputdict),
                    0,
                    0,
                    0,
                ],
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            wfndosdict['system']['noncolin'] = True
            wfndosdict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(wfndosdict, jmespath.search('dos.args', self.inputdict))

        return QeGrammar().write(wfndosdict)

    @property
    def dos(self):
        dosdict: dict = {
            'dos': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'fildos': "'dos.dat'",
            }
        }

        # Update if needed. 
        update_dict(dosdict, jmespath.search('dos.args', self.inputdict))

        file_string: str =  NamelistGrammar().write(dosdict)

        return file_string

    @property
    def job_wfndos(self):
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'dos.job_wfndos_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

rm -rf ./tmp
cp -r ../scf/tmp ./tmp
{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < wfn.in &> wfn.in.out
'''
        return file_string

    @property
    def job_dos(self):
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'dos.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}dos.x -pd .true. {scheduler.get_exec_infix()} < dos.in &> dos.in.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './dos/wfn.in': self.wfndos,
            './dos/dos.in': self.dos,
            './dos/job_wfn.sh': self.job_wfndos,
            './dos/job_dos.sh': self.job_dos,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './dos/job_wfn.sh',
            './dos/job_dos.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './dos',
        ]
    
    def plot(self, **kwargs):
            DosPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion