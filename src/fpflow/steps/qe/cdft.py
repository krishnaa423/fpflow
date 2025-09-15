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
import numpy as np 
#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class QeCdftStep(Step):
    @property
    def cdft(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = int(qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))
        cond_band: int = jmespath.search('cdft.excited_cond_band', self.inputdict)
        val_band: int = jmespath.search('cdft.removed_val_band', self.inputdict)
        

        qedict: dict = {
            'control': {
                'outdir': './tmp',
                'prefix': 'struct',
                'pseudo_dir': './pseudos/qe',
                'calculation': jmespath.search('cdft.type', self.inputdict),
                'tprnfor': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict),
                'nbnd': int(max_val_bands + cond_band),
                'occupations': 'from_input',
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
            },
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            qedict['system']['noncolin'] = True
            qedict['system']['lspinorb'] = True

            # Set the occupations array. 
            occupations = np.zeros((max_val_bands + cond_band), dtype='i4')
            occupations[:max_val_bands] = 1       
            occupations[max_val_bands + val_band] = 0       # remove val band state. 
            occupations[max_val_bands + (cond_band-1)] = 1  # Populate cond band state. 
            qedict['occupations'] = occupations.tolist()
        else:
            # Set the occupations array. No spin-orbit case. 
            occupations = np.zeros((max_val_bands + cond_band), dtype='i4')
            occupations[:max_val_bands] = 2      
            occupations[max_val_bands + val_band] = 0       # remove val band state. 
            occupations[max_val_bands + (cond_band-1)] = 2  # Populate cond band state. 
            qedict['occupations'] = occupations.tolist()

        # Update if needed. 
        update_dict(qedict, jmespath.search('cdft.args', self.inputdict))

        return QeGrammar().write(qedict)

    @property
    def job_cdft(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'cdft.job_info')

        save_final_cell_parameters_str = "awk '/Begin final coordinates/ {end_flag=1; next} end_flag && /CELL_PARAMETERS/ {cell_flag=1; next} /End final coordinates/ {end_flag=0} end_flag && cell_flag {print; if (length==0) cell_flag=0 }' cdft.in.out > relaxed_cell_parameters.txt"
        save_final_atomic_positions_str = "awk '/Begin final coordinates/ {end_flag=1; next} end_flag && /ATOMIC_POSITIONS/ {pos_flag=1; next} /End final coordinates/ {end_flag=0}  end_flag && pos_flag { print $1, $2, $3, $4 }' cdft.in.out > relaxed_atomic_positions.txt"
        update_coord: bool = jmespath.search('cdft.update_coord', self.inputdict)


        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < cdft.in &> cdft.in.out

cp ./tmp/struct.save/data-file-schema.xml ./cdft.xml

# Copy the end atomic positions and cell parameters. 
{save_final_cell_parameters_str if update_coord else ""}
{save_final_atomic_positions_str if update_coord else ""}

# Update from relax.
{"fpflow generator --create" if update_coord else ""}
'''
        
        return file_string
    
    @property
    def file_contents(self) -> dict:
        return {
            'cdft.in': self.cdft,
            'job_cdft.sh': self.job_cdft,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return ['./job_cdft.sh']

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './cdft.in',
            './job_cdft.sh',
            './tmp',
            './cdft.out',
            './cdft.in.out',
            './cdft.xml',
            './relaxed_cell_parameters.txt',
            './relaxed_atomic_positions.txt',
        ]

#endregion