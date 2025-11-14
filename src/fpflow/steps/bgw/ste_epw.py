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
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.io.logging import get_logger
from fpflow.plots.ste_epw import EpwStePlot

#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class EpwSteStep(Step):
    @property
    def bands_skipped_string(self) -> str:
        #TODO: Copied this code, but need to refactor to make it simple. 

        bands_skipped = None 

        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = int(qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))

        # Populate list.
        if bands_skipped is None:
            bands_skipped = []
            epw_val_bands = jmespath.search('elph.val_bands', self.inputdict)
            total_val_bands = max_val_bands
            epw_cond_bands = jmespath.search('elph.cond_bands', self.inputdict)
            wfn_cond = jmespath.search('cond_bands', self.wfn_options)

            if epw_val_bands!= total_val_bands:
                temp = (1, total_val_bands - epw_val_bands)
                bands_skipped.append(temp)

            if epw_cond_bands!= wfn_cond and epw_cond_bands < wfn_cond:
                temp = (total_val_bands+ epw_cond_bands + 1, wfn_cond + total_val_bands)
                bands_skipped.append(temp)

            if len(bands_skipped)==0:
                bands_skipped = None

        # Populate string. 
        bands_skipped_str = ''
        exclude_bands_str = None
        if bands_skipped is not None:
            num_bands_skipped = len(bands_skipped)
            exclude_bands_str = "'exclude_bands="
            
            for bands_idx, bands in enumerate(bands_skipped):
                exclude_bands_str += f'{bands[0]}:{bands[1]}'
                if bands_idx!=num_bands_skipped-1: exclude_bands_str += ','
                
            exclude_bands_str += "'"
            
            bands_skipped_str = 'bands_skipped=' + exclude_bands_str
        
        return exclude_bands_str

    @property
    def ste_epw(self) -> str:
        # wfnlink. 
        wfn_link: str = jmespath.search('elph.nscf_link', self.inputdict)
        wfn_names: list[str] = jmespath.search(f'nscf.list[*].name', self.inputdict)
        wfn_index: int = wfn_names.index(wfn_link)
        self.wfn_options: dict = jmespath.search(f'nscf.list[{wfn_index}]', self.inputdict)

        epwdict: dict = {
            'inputepw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'dvscf_dir': "'./save'",

                'nk1': jmespath.search('ste.kgrid[0]', self.inputdict),
                'nk2': jmespath.search('ste.kgrid[1]', self.inputdict),
                'nk3': jmespath.search('ste.kgrid[2]', self.inputdict),
                'nq1': jmespath.search('ste.qgrid[0]', self.inputdict),
                'nq2': jmespath.search('ste.qgrid[1]', self.inputdict),
                'nq3': jmespath.search('ste.qgrid[2]', self.inputdict),
                'nkf1': jmespath.search('ste.kgrid[0]', self.inputdict),
                'nkf2': jmespath.search('ste.kgrid[1]', self.inputdict),
                'nkf3': jmespath.search('ste.kgrid[2]', self.inputdict),
                'nqf1': jmespath.search('ste.qgrid[0]', self.inputdict),
                'nqf2': jmespath.search('ste.qgrid[1]', self.inputdict),
                'nqf3': jmespath.search('ste.qgrid[2]', self.inputdict),
                'nbndsub': jmespath.search('elph.val_bands', self.inputdict) + jmespath.search('elph.cond_bands', self.inputdict),
                
                'elph': '.true.',
                # 'epbwrite': '.true.',
                # 'epbread': '.false.',
                'epwwrite': '.false.',
                'epwread': '.true.',
                'lpolar': '.true.',

                'exciton': '.true.',
                'explrn': '.true.',
                'negnv_explrn': jmespath.search('ste.nxct', self.inputdict),
                'nbndv_explrn': jmespath.search('ste.val_bands', self.inputdict),
                'nbndc_explrn': jmespath.search('ste.cond_bands', self.inputdict) - 1, # -1 for the fine grid thing. 
                'ethrdg_plrn': jmespath.search('ste.max_error', self.inputdict),
                'init_plrn': 5,
                'niter_plrn': jmespath.search('ste.max_steps', self.inputdict),

                'wannierize': '.false.',
                
                'iprint': 2,
            }
        }

        # Add bands_skipped.
        a = self.bands_skipped_string
        if self.bands_skipped_string is not None and self.bands_skipped_string!='' : 
            epwdict['inputepw']['bands_skipped'] = self.bands_skipped_string

        # Update if needed. 
        update_dict(epwdict, jmespath.search('ste.args', self.inputdict))

        # If projections provided, delete auto_projections and scdm_proj.
        if 'proj(1)' in epwdict['inputepw'].keys() and epwdict['inputepw']['wannierize'] == '.true.': 
            del epwdict['inputepw']['auto_projections']
            del epwdict['inputepw']['scdm_proj']

        return NamelistGrammar().write(epwdict)
    
    @property
    def ste_epw_plrn_plot(self) -> str:
        # wfnlink. 
        wfn_link: str = jmespath.search('elph.nscf_link', self.inputdict)
        wfn_names: list[str] = jmespath.search(f'nscf.list[*].name', self.inputdict)
        wfn_index: int = wfn_names.index(wfn_link)
        self.wfn_options: dict = jmespath.search(f'nscf.list[{wfn_index}]', self.inputdict)

        epwdict: dict = {
            'inputepw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'dvscf_dir': "'./save'",

                'nk1': jmespath.search('ste.kgrid[0]', self.inputdict),
                'nk2': jmespath.search('ste.kgrid[1]', self.inputdict),
                'nk3': jmespath.search('ste.kgrid[2]', self.inputdict),
                'nq1': jmespath.search('ste.qgrid[0]', self.inputdict),
                'nq2': jmespath.search('ste.qgrid[1]', self.inputdict),
                'nq3': jmespath.search('ste.qgrid[2]', self.inputdict),
                'nkf1': jmespath.search('ste.kgrid[0]', self.inputdict),
                'nkf2': jmespath.search('ste.kgrid[1]', self.inputdict),
                'nkf3': jmespath.search('ste.kgrid[2]', self.inputdict),
                'nqf1': jmespath.search('ste.qgrid[0]', self.inputdict),
                'nqf2': jmespath.search('ste.qgrid[1]', self.inputdict),
                'nqf3': jmespath.search('ste.qgrid[2]', self.inputdict),
                'nbndsub': jmespath.search('elph.val_bands', self.inputdict) + jmespath.search('elph.cond_bands', self.inputdict),
                
                'elph': '.true.',
                # 'epbwrite': '.true.',
                # 'epbread': '.false.',
                'epwwrite': '.false.',
                'epwread': '.true.',
                'lpolar': '.true.',

                'exciton': '.true.',
                'explrn': '.true.',
                'negnv_explrn': jmespath.search('ste.nxct', self.inputdict),
                'nbndv_explrn': jmespath.search('ste.val_bands', self.inputdict),
                'nbndc_explrn': jmespath.search('ste.cond_bands', self.inputdict)- 1, # -1 for the fine grid thing.
                'ethrdg_plrn': jmespath.search('ste.max_error', self.inputdict),
                'init_plrn': 5,
                'niter_plrn': jmespath.search('ste.max_steps', self.inputdict),
                'plot_explrn_e': '.true.',
                'plot_explrn_h': '.true.',

                'wannierize': '.false.',
                
                'iprint': 2,
            }
        }

        # Add bands_skipped.
        a = self.bands_skipped_string
        if self.bands_skipped_string is not None and self.bands_skipped_string!='' : 
            epwdict['inputepw']['bands_skipped'] = self.bands_skipped_string

        # Update if needed. 
        update_dict(epwdict, jmespath.search('ste.args', self.inputdict))

        # If projections provided, delete auto_projections and scdm_proj.
        if 'proj(1)' in epwdict['inputepw'].keys() and epwdict['inputepw']['wannierize'] == '.true.': 
            del epwdict['inputepw']['auto_projections']
            del epwdict['inputepw']['scdm_proj']

        return NamelistGrammar().write(epwdict)

    @property
    def job_ste_epw(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'ste.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

# Prereq: xctph_epw run must be done to generate the necessary files and folders. 
{scheduler.get_exec_prefix()}epw.x -npool {scheduler.nk} < ste.in  &> ste.in.out 
'''
        return file_string
    
    @property
    def job_ste_epw_plrn_plot(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'ste.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

# Prereq: xctph_epw run must be done to generate the necessary files and folders.
{scheduler.get_exec_prefix()}epw.x -npool {scheduler.nk} < ste_wf_plot.in  &> ste_wf_plot.in.out 
'''
        return file_string


    @property
    def file_contents(self) -> dict:
        return {
            './zd_epw/ste.in': self.ste_epw,
            './zd_epw/ste_wf_plot.in': self.ste_epw_plrn_plot,
            './zd_epw/job_ste.sh': self.job_ste_epw,
            './zd_epw/job_ste_wf_plot.sh': self.job_ste_epw_plrn_plot,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './zd_epw/job_ste.sh',
            './zd_epw/job_ste_wf_plot.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './zd_epw',
        ]
    
    def plot(self, **kwargs):
        EpwStePlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion