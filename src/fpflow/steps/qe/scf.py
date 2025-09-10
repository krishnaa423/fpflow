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
                'pseudo_dir': './pseudos/qe',
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

{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < scf.in &> scf.in.out

cp ./tmp/struct.save/data-file-schema.xml ./scf.xml
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'scf.in': self.scf,
            'job_scf.sh': self.job_scf,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return ['./job_scf.sh']

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './scf.in',
            './job_scf.sh',
            './tmp',
            './scf.in.out',
            './scf.xml',
            './pseudos/qe',
            './atoms_*.xsf',        # Removes structure files created by Struct object. 
        ]
#endregion