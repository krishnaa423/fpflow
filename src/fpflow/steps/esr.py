#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.nestedstep_base import NestedBaseStep
import jmespath 
from fpflow.schedulers.scheduler import Scheduler
import copy 
import numpy as np 
from glom import Path 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class EsrStep(NestedBaseStep):
    def __init__(self, **kwargs):
        super().__init__(steptag='esr', **kwargs)

    def get_band_changes(self):
        supercell_size: np.ndarray = np.array(jmespath.search('esr.supercell_size', self.inputdict))
        multiplier: int = int(np.prod(supercell_size).item())

        output: list = [
            {
                'path': 'dftelbands.cond_bands',
                'value': jmespath.search('dftelbands.cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'wfn.cond_bands',
                'value': jmespath.search('wfn.cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'wfnq.cond_bands',
                'value': jmespath.search('wfnq.cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'wfn.parabands_cond_bands',
                'value': jmespath.search('wfn.parabands_cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'gw.epsilon.cond_bands',
                'value': jmespath.search('gw.epsilon.cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'gw.sigma.cond_bands',
                'value': jmespath.search('gw.sigma.cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'gw.sigma.conv_cond_bands',
                'value': jmespath.search('gw.sigma.conv_cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'gw.sigma.val_bands',
                'value': jmespath.search('gw.sigma.val_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'gw.gwelbands.val_bands',
                'value': jmespath.search('gw.gwelbands.val_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'bse.absorption.val_bands',
                'value': jmespath.search('bse.absorption.val_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'bse.absorption.cond_bands',
                'value': jmespath.search('bse.absorption.cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'bse.absorption.nxct',
                'value': jmespath.search('bse.absorption.nxct', self.inputdict)*multiplier,  
            },
            # Now the common.
            {
                'path': 'common.dft.val_bands',
                'value': jmespath.search('common.dft.val_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'common.dft.cond_bands',
                'value': jmespath.search('common.dft.cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'common.gw.val_bands',
                'value': jmespath.search('common.gw.val_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'common.gw.plot_val_bands',
                'value': jmespath.search('common.gw.plot_val_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'common.gw.cond_bands',
                'value': jmespath.search('common.gw.cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'common.gw.para_cond_bands',
                'value': jmespath.search('common.gw.para_cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'common.gw.conv_cond_bands',
                'value': jmespath.search('common.gw.conv_cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'common.bse.val_bands',
                'value': jmespath.search('common.bse.val_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'common.bse.cond_bands',
                'value': jmespath.search('common.bse.cond_bands', self.inputdict)*multiplier,  
            },
            {
                'path': 'common.bse.nxct',
                'value': jmespath.search('common.bse.nxct', self.inputdict)*multiplier,  
            },
        ]

        return output 

    def get_kgrid_changes(self):
        gamma_kgrid: list = [1, 1, 1]

        output: list = [
            {
                'path': 'scf.kgrid',
                'value': gamma_kgrid,
            },
            {
                'path': 'scf.args.system.nosym',
                'value': True,
            },
            {
                'path': 'dfpt.qgrid',
                'value': gamma_kgrid,
            },
            {
                'path': 'epw.job_info.nk',
                'value': 1,
            },
            {
                'path': 'wfn.kgrid',
                'value': gamma_kgrid,
            },
            {
                'path': 'wfn.sym',
                'value': False,
            },
            {
                'path': 'wfnq.kgrid',
                'value': gamma_kgrid,
            },
            {
                'path': 'wfnq.sym',
                'value': False,
            },
            {
                'path': 'bseq.qgrid',
                'value': gamma_kgrid,
            },
            {
                'path': 'bseq.sym',
                'value': False,
            },
            {
                'path': 'bse.plotxct.supercell_size',
                'value': gamma_kgrid,
            },
            ## Now the common. 
            {
                'path': 'common.sym',
                'value': False,
            },
            {
                'path': 'common.kgrid',
                'value': gamma_kgrid,
            },
            {
                'path': 'common.nkpts',
                'value': 1,
            },
            {
                'path': 'common.qgrid',
                'value': gamma_kgrid,
            },
            {
                'path': 'common.epw_nodes',
                'value': 1,
            },
            {
                'path': 'common.bse.supercell_size',
                'value': gamma_kgrid,
            },
        ]

        return output

    @property
    def folder_inputdict_changes_map(self) -> dict:
        '''
        Return:
            {
                'folder1': [
                    {
                        'path': 'some.path1',
                        'value': 12,
                    },
                    {
                        'path': 'some.path2',
                        'value': 45,
                    },
                ],
                'folder2': [
                    {
                        'path': 'some.path1',
                        'value': 123,
                    },
                    {
                        'path': 'some.path2',
                        'value': 345,
                    },
                ],
            }
        '''
        max_steps = jmespath.search('esr.max_steps', self.inputdict)

        output = {}

        generator_pre_steps: list = ['pseudos_qe']
        
        generator_non_initial_steps: list = [
            'scf_qe',
            'dfpt_qe',
            'wfn_qe',
            'epw_qe',
            'wfnq_qe',
            'epsilon_bgw',
            'sigma_bgw',
            'bseq_bgw',
            'esf',
        ]
        generator_initial_steps: list = copy.deepcopy(generator_non_initial_steps)
        generator_initial_steps.insert(0, 'cdft_qe')
        generator_non_initial_steps.insert(0, 'esr_gen')
        generator_post_steps: list = ['create_script', 'run_script', 'remove_script', 'plot_script', 'interactive_script']
        manager_non_initial_steps: list = [
            'scf_qe',
            'dfpt_qe',
            'wfn_qe',
            'epw_qe',
            'wfnq_qe',
            'epsilon_bgw',
            'sigma_bgw',
            'bseq_bgw',
            'esf',
        ]
        manager_initial_steps: list = copy.deepcopy(manager_non_initial_steps)
        manager_initial_steps.insert(0, 'cdft_qe')
        manager_non_initial_steps.insert(0, 'esr_gen')

        kgrid_changes_list: list = self.get_kgrid_changes()
        band_changes_list: list = self.get_band_changes()

        for step_idx in range(max_steps):
            generator_steps: list = generator_initial_steps if step_idx==0 else generator_non_initial_steps
            manager_steps: list = manager_initial_steps if step_idx==0 else manager_non_initial_steps
            
            output[f'./esr/dset_{step_idx}/'] = [
                {
                    'path': 'structures.active_idx',
                    'value': 0,
                },
                {
                    'path': 'structures.list[0].supercell_size',
                    'value': jmespath.search('esr.supercell_size', self.inputdict),
                },
                {
                    'path': 'generator.pre_steps',
                    'value': generator_pre_steps,
                },
                {
                    'path': 'generator.steps',
                    'value': generator_steps,
                },
                {
                    'path': 'generator.post_steps',
                    'value': generator_post_steps,
                },
                {
                    'path': 'manager.steps',
                    'value': manager_steps,
                },
                {
                    'path': 'generator.dest_dir',
                    'value': './',
                }
            ] + kgrid_changes_list + band_changes_list

        return output

    @property
    def extra_filecontents(self) -> dict:
        output: dict = super().extra_filecontents

        return output
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_esr.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './esr',
            './job_esr.sh',
            './script_esr.py',
            './script_esr.py.out',
        ]

#endregion