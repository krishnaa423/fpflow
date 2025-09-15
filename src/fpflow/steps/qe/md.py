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
class QeMdStep(Step):
    @property
    def md(self):
        qestruct = QeStruct.from_inputdict(self.inputdict)
        qestruct.struct_idx = jmespath.search('md.structure_index', self.inputdict)

        qedict: dict = {
            'control': {
                'outdir': './tmp',
                'prefix': 'struct',
                'pseudo_dir': './pseudos/qe',
                'calculation': jmespath.search('md.type', self.inputdict),
                'tprnfor': True,
                'dt': jmespath.search('md.time_step', self.inputdict),
                'nstep': jmespath.search('md.nsteps', self.inputdict),
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict),
                'nosym': True,
            },
            'electrons': {},
            'ions': {
                'pot_extrapolation': 'second-order',
                'wfc_extrapolation': 'second-order',
                'ion_temperature': 'initial',
                'tempw': jmespath.search('md.temp', self.inputdict),
            },
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'automatic',
                'data': [
                    jmespath.search('md.kgrid[0]', self.inputdict),
                    jmespath.search('md.kgrid[1]', self.inputdict),
                    jmespath.search('md.kgrid[2]', self.inputdict),
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
        update_dict(qedict, jmespath.search('md.args', self.inputdict))

        return QeGrammar().write(qedict)

    @property
    def job_md(self):
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'md.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < md.in &> md.in.out

cp ./tmp/struct.save/data-file-schema.xml ./md.xml
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'md.in': self.md,
            'job_md.sh': self.job_md,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return ['./job_md.sh']

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './md.in',
            './job_md.sh',
            './tmp',
            './md.in.out',
            './md.xml',
            './pseudos/qe',
            './atoms_*.xsf',        # Removes structure files created by Struct object. 
        ]
#endregion