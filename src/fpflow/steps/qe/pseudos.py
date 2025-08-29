#region modules
from fpflow.steps.step import Step
from typing import List 
from fpflow.inputs.grammars.qe import QeGrammar
from fpflow.structure.qe.qe_struct import QeStruct
import jmespath
import os 
from fpflow.inputs.inputyaml import InputYaml
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QePseudosStep(Step):
    def generate_pseudos(self):
        qestruct = QeStruct.from_inputdict(self.inputdict)
        xc = jmespath.search('scf.xc', self.inputdict)
        is_soc = jmespath.search('scf.is_spinorbit', self.inputdict)

        paths, filenames = qestruct.get_pseudos_list(xc=xc, is_soc=is_soc)

        os.system('mkdir -p pseudos/qe')
        for path, filename in zip(paths, filenames):
            os.system(f'cp {path} ./pseudos/qe/{filename}')

    @property
    def file_contents(self) -> dict:
        return {
            'script_pseudos.py': f'''#!/usr/bin/env python3

from fpflow.inputs.inputyaml import InputYaml
from fpflow.steps.qe.pseudos import QePseudosStep

inputdict: dict = InputYaml.from_yaml_file('./input.yaml').inputdict
QePseudosStep(inputdict).generate_pseudos()
''',
            'job_pseudos.sh': f'''#!/bin/bash

python3 script_pseudos.py
'''
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return ['./job_pseudos.sh']

    @property
    def save_inodes(self) -> List[str]:
        return ['./pseudos']
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './pseudos',
            './script_pseudos.py',
            './job_pseudos.sh',
        ]
#endregion