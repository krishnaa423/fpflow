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
from fpflow.plots.xctph_epw import EpwXctphPlot

#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class EpwXctphStep(Step):
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

        # wfnlink. 
        wfn_link: str = jmespath.search('elph.nscf_link', self.inputdict)
        wfn_names: list[str] = jmespath.search(f'nscf.list[*].name', self.inputdict)
        wfn_index: int = wfn_names.index(wfn_link)
        self.wfn_options: dict = jmespath.search(f'nscf.list[{wfn_index}]', self.inputdict)

        # Populate list.
        if bands_skipped is None:
            bands_skipped = []
            epw_val_bands = jmespath.search('xctph.val_bands', self.inputdict)
            total_val_bands = max_val_bands
            epw_cond_bands = jmespath.search('xctph.cond_bands', self.inputdict)
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
    def xctph_epw(self) -> str:
        epwdict: dict = {
            'inputepw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'dvscf_dir': "'./save'",
                
                'nk1': jmespath.search('xctph.kgrid[0]', self.inputdict),
                'nk2': jmespath.search('xctph.kgrid[1]', self.inputdict),
                'nk3': jmespath.search('xctph.kgrid[2]', self.inputdict),
                'nq1': jmespath.search('xctph.qgrid[0]', self.inputdict),
                'nq2': jmespath.search('xctph.qgrid[1]', self.inputdict),
                'nq3': jmespath.search('xctph.qgrid[2]', self.inputdict),
                'nkf1': jmespath.search('xctph.kgrid[0]', self.inputdict),
                'nkf2': jmespath.search('xctph.kgrid[1]', self.inputdict),
                'nkf3': jmespath.search('xctph.kgrid[2]', self.inputdict),
                'nqf1': jmespath.search('xctph.qgrid[0]', self.inputdict),
                'nqf2': jmespath.search('xctph.qgrid[1]', self.inputdict),
                'nqf3': jmespath.search('xctph.qgrid[2]', self.inputdict),
                'nbndsub': jmespath.search('xctph.val_bands', self.inputdict) + jmespath.search('xctph.cond_bands', self.inputdict),
                
                'elph': '.true.',
                # 'epbwrite': '.true.',
                # 'epbread': '.false.',
                'epwwrite': '.true.',
                'epwread': '.false.',
                'lpolar': '.true.',

                'exciton': '.true.',
                'explrn': '.false.',
                'negnv_explrn': jmespath.search('xctph.nxct', self.inputdict),
                'nbndv_explrn': jmespath.search('xctph.val_bands', self.inputdict),
                'nbndc_explrn': jmespath.search('xctph.cond_bands', self.inputdict),
                
                'wannierize': '.true.',
                'wannier_plot':'.true.',
                'wannier_plot_supercell': f"{jmespath.search('xctph.qgrid[0]', self.inputdict)*2} {jmespath.search('xctph.qgrid[1]', self.inputdict)*2} {jmespath.search('xctph.qgrid[2]', self.inputdict)*2}",
                'auto_projections': '.true.',
                'scdm_proj': '.true.',
                
                'iprint': 2,
            }
        }

        # Add bands_skipped.
        a = self.bands_skipped_string
        if self.bands_skipped_string is not None and self.bands_skipped_string!='' : 
            epwdict['inputepw']['bands_skipped'] = self.bands_skipped_string

        # Update if needed. 
        update_dict(epwdict, jmespath.search('xctph.args', self.inputdict))

        return NamelistGrammar().write(epwdict)

    @property
    def job_xctph_epw(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'xctph.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

rm -rf ./tmp
ln -sf ../{self.wfn_options["name"]}/tmp ./tmp
ln -sf ../dfpt/save ./save
ln -sf ../bseq ./eigv
{scheduler.get_exec_prefix()}epw.x {scheduler.get_exec_infix()} < xctph.in  &> xctph.in.out 
'''
        return file_string


    @property
    def file_contents(self) -> dict:
        return {
            './xctph_epw/xctph.in': self.xctph_epw,
            './xctph_epw/job_xctph.sh': self.job_xctph_epw,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './xctph_epw/job_xctph.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './xctph_epw',
        ]

    def plot(self, **kwargs):
        EpwXctphPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion