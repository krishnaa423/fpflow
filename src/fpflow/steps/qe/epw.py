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
#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class QeEpwStep(Step):
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
            abs_val_bands = jmespath.search('bse.absorption.val_bands', self.inputdict)
            total_val_bands = max_val_bands
            abs_cond_bands = jmespath.search('bse.absorption.cond_bands', self.inputdict)
            wfn_cond = jmespath.search('wfn.cond_bands', self.inputdict)

            if abs_val_bands!= total_val_bands:
                temp = (1, total_val_bands - abs_val_bands)
                bands_skipped.append(temp)

            if abs_cond_bands!= wfn_cond and abs_cond_bands < wfn_cond:
                temp = (total_val_bands+ abs_cond_bands + 1, wfn_cond + total_val_bands)
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
    def epw(self) -> str:
        epwdict: dict = {
            'inputepw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'nk1': jmespath.search('dfpt.qgrid[0]', self.inputdict),
                'nk2': jmespath.search('dfpt.qgrid[1]', self.inputdict),
                'nk3': jmespath.search('dfpt.qgrid[2]', self.inputdict),
                'nq1': jmespath.search('dfpt.qgrid[0]', self.inputdict),
                'nq2': jmespath.search('dfpt.qgrid[1]', self.inputdict),
                'nq3': jmespath.search('dfpt.qgrid[2]', self.inputdict),
                'nkf1': jmespath.search('dfpt.qgrid[0]', self.inputdict),
                'nkf2': jmespath.search('dfpt.qgrid[1]', self.inputdict),
                'nkf3': jmespath.search('dfpt.qgrid[2]', self.inputdict),
                'nqf1': jmespath.search('dfpt.qgrid[0]', self.inputdict),
                'nqf2': jmespath.search('dfpt.qgrid[1]', self.inputdict),
                'nqf3': jmespath.search('dfpt.qgrid[2]', self.inputdict),
                'nbndsub': jmespath.search('bse.absorption.val_bands', self.inputdict) + jmespath.search('bse.absorption.cond_bands', self.inputdict),
                'dvscf_dir': "'./save'",
                'elph': '.true.',
                'epbwrite': '.true.',
                'epbread': '.false.',
                'prtgkk': '.false.',
                'wannierize': '.true.',
                'auto_projections': '.true.',
                'scdm_proj': '.true.',
                # 'temps': '300.0',
                # 'verbosity': 1,
            }
        }

        # Add bands_skipped.
        a = self.bands_skipped_string
        if self.bands_skipped_string is not None and self.bands_skipped_string!='' : 
            epwdict['inputepw']['bands_skipped'] = self.bands_skipped_string

        # Update if needed. 
        update_dict(epwdict, jmespath.search('epw.args', self.inputdict))

        return NamelistGrammar().write(epwdict)

    @property
    def job_epw(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'epw.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}epw.x {scheduler.get_exec_infix()} < epw.in  &> epw.in.out 
cp ./wfn.xml ./save/wfn.xml
cp ./tmp/*epb* ./save/
'''
        return file_string


    @property
    def file_contents(self) -> dict:
        return {
            'epw.in': self.epw,
            'job_epw.sh': self.job_epw,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_epw.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './epw.in',
            './job_epw.sh',
            './tmp',
            './struct*',
            './decay*',
            './EPW.bib',
            './epwdata.fmt',
            './selecq.fmt',
            './vmedata.fmt',
            './epw.in.out',
            './epw.out',
            './crystal.fmt',
        ]
#endregion