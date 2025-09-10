#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.inputs.grammars.namelist import NamelistGrammar
import jmespath
from fpflow.io.update import update_dict
from fpflow.io.logging import get_logger
from fpflow.schedulers.scheduler import Scheduler
from importlib.util import find_spec
import glob 
from ase.dft.kpoints import get_special_points
from fpflow.structure.struct import Struct
from ase import Atoms
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.inputs.grammars.qe import QeGrammar
from fpflow.plots.phbands import PhonopyPlot
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QePhonopyStep(Step):
    @property
    def phonopy_scf_prefix(self) -> str:
        qestruct = QeStruct.from_inputdict(self.inputdict)

        qeprefixdict: dict = {
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
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            qeprefixdict['system']['noncolin'] = True
            qeprefixdict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(qeprefixdict, jmespath.search('scf.args', self.inputdict))

        return QeGrammar().write(qeprefixdict)
    
    @property
    def phonopy_scf_suffix(self) -> str:
        qesuffixdict: dict = {
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
            qesuffixdict['system']['noncolin'] = True
            qesuffixdict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(qesuffixdict, jmespath.search('scf.args', self.inputdict))

        return QeGrammar().write(qesuffixdict)
    
    def get_phonopy_bandpath(self) -> str:
        #TODO: Copy pasted. Can refactor this. 
        struct: Struct = Struct.from_inputdict(self.inputdict)
        atoms: Atoms = struct.atoms[struct.struct_idx]
        sc_map = get_special_points(atoms.cell)
        sc_labels = jmespath.search('kpath.special_points[*]', self.inputdict)
        npoints_segment = jmespath.search('kpath.npoints_segment', self.inputdict)
        output = ''

        output += '--band=" '

        for sc_label in sc_labels:
            for col in sc_map[sc_label]:
                output += f' {col:15.10f} '

        output += ' " '

        output += f' --band-points={npoints_segment} '

        return output

    @property
    def job_phonopy(self) -> str:
        #TODO: Copy pasted. Can refactor this. 

        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'dfpt.job_info')

        qdim: list = jmespath.search('dfpt.qgrid[*]', self.inputdict)
        str_2_f(self.phonopy_scf_prefix, 'phonopy_scf_prefix')
        str_2_f(self.phonopy_scf_suffix, 'phonopy_scf_suffix')
        os.system(f'phonopy --qe -d --dim="{qdim[0]} {qdim[1]} {qdim[2]}" -c scf.in')
        sc_files = glob.glob('supercell-*')
        for sc_file in sc_files:
            os.system(f'cat phonopy_scf_prefix {sc_file} phonopy_scf_suffix >| phonopy-{sc_file}')


        # Create supercell job.
        phonopy_bandpath_str: str = self.get_phonopy_bandpath()
        files = glob.glob('phonopy-supercell-*')
        start_idx = 0
        stop_idx = len(files)
        debug_str: str = '\n'
        files_bashvar_str: str = '\nfiles=('
        files_args_str: str = ''
        for file_idx, file in enumerate(files): 
            files_bashvar_str += f'"{file}" '
            files_args_str += f' {file}.out '
            debug_str += f'#idx: {file_idx}, filename: {file}\n'
        files_bashvar_str += ')\n\n'
        debug_str += '\n\n'
        file_variable = '${files[$i]}'

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{debug_str}

{files_bashvar_str}

start={start_idx}
stop={stop_idx-1}

for (( i=$start; i<=$stop; i++ )); do
{scheduler.get_exec_prefix()}pw.x < {file_variable} &> {file_variable}.out
done

# Post processing. This should create FORCE_SETS and phonon bands. 
phonopy -f {files_args_str}
phonopy --qe -c scf.in {phonopy_bandpath_str} --dim="{qdim[0]} {qdim[1]} {qdim[2]}"
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'phonopy_scf_prefix': self.phonopy_scf_prefix,
            'phonopy_scf_suffix': self.phonopy_scf_suffix,
            'job_phonopy.sh': self.job_phonopy,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_phonopy.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './phonopy_scf_prefix',
            './phonopy_disp.yaml',
            './phonopy_scf_suffix',
            './phonopy-supercell*',
            './supercell*',
            './job_phonopy.sh',
            './FORCE_SETS',
            './phonopy.yaml',
            './band.yaml',
        ]
    
    def plot(self):
        PhonopyPlot().save_plot()

#endregion