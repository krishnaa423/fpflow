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
from fpflow.structure.kpts import Kpts
import math 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeDfptStep(Step):
    @property
    def create_save(self) -> str:
        file_string = ''

        pkg_dir = os.path.dirname(find_spec('fpflow').origin)
        path = os.path.join(pkg_dir, 'data', 'create_save.txt')
        with open(path, 'r') as f:
            file_string = f.read()

        return file_string

    @property
    def dfpt(self) -> str:
        dfptdict: dict = {
            'inputph': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'ldisp': '.true.',
                'nq1': jmespath.search('dfpt.qgrid[0]', self.inputdict),
                'nq2': jmespath.search('dfpt.qgrid[1]', self.inputdict),
                'nq3': jmespath.search('dfpt.qgrid[2]', self.inputdict),
                'fildyn': "'struct.dyn'",
                'tr2_ph': jmespath.search('dfpt.conv_thr', self.inputdict),
                'fildvscf': "'dvscf'",
            }
        }

        # Update if needed. 
        update_dict(dfptdict, jmespath.search('dfpt.args', self.inputdict))

        return NamelistGrammar().write(dfptdict)
    
    @property
    def job_dfpt(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'dfpt.job_info')
        auto_assign_ni_and_nodes: bool = jmespath.search('dfpt.auto_assign_ni_and_nodes', self.inputdict)

        if auto_assign_ni_and_nodes is not None and auto_assign_ni_and_nodes:
            # Assign ni. 
            nqirr: int = len(Kpts.from_kgrid(
                kgrid = jmespath.search('dfpt.qgrid', self.inputdict),
                is_reduced=True,
            ).kpts)
            scheduler.ni = nqirr

            # Assign nodes. 
            cores_per_node: int = jmespath.search('cores', scheduler.node_info)
            requested_nodes = scheduler.nodes
            requested_tasks = scheduler.ntasks
            assigned_nodes: int = 1
            assigned_tasks: int = 1
            #If tasks are ge then qpts. 
            if requested_tasks >= nqirr:
                assert requested_tasks <= requested_nodes*cores_per_node, 'requested tasks is greater than total cores.'

                assigned_tasks = (requested_tasks // nqirr ) * nqirr
            #If tasks are lt than kpts. Atleast set number of tasks as nqirr. 
            else:
                assigned_tasks = nqirr
            
            assigned_nodes = math.ceil(requested_tasks/cores_per_node)

            # Finally update scheduler.
            scheduler.nodes = assigned_nodes
            scheduler.ntasks = assigned_tasks
            scheduler.ni = nqirr

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}ph.x {scheduler.get_exec_infix()} < dfpt.in &> dfpt.in.out

python3 ./create_save.py
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            'dfpt.in': self.dfpt,
            'job_dfpt.sh': self.job_dfpt,
            'create_save.py': self.create_save,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return ['./job_dfpt.sh']

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './dfpt*.in',
            './dfpt*.in.out',
            './create_save.py',
            './job_dfpt.sh',
            './tmp',
            './out*',
            './save',
            './struct.dyn*',
        ]

#endregion