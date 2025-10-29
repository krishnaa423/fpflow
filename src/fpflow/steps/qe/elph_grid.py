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
from fpflow.structure.kpath import Kpath
from fpflow.inputs.grammars.wannier import WannierGrammar
from fpflow.plots.elph_grid import EpwElphGridPlot

#endregion

#region variables
logger = get_logger()
#endregion

#region functions
#endregion

#region classes
class EpwElphGridStep(Step):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Set options.
        wfnlink: str = jmespath.search('elph.nscf_link', self.inputdict)
        wfn_names: list[str] = jmespath.search(f'nscf.list[*].name', self.inputdict)
        wfn_index: int = wfn_names.index(wfnlink)
        self.wfn_options: dict = jmespath.search(f'nscf.list[{wfn_index}]', self.inputdict)

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

    def add_wannier_data(self, epwdict: dict):
        # Kpath.
        kpath_data: list = Kpath.from_yamlfile().wannierpath_list
        
        wannierdict: dict = {
            'bands_plot': True,
            'bands_num_points': jmespath.search('kpath.npoints_segment', self.inputdict),
            'write_hr': True,
            'use_ws_distance': True,
            'write_xyz': True,
            'write_u_matrices': True,
            'write_tb': True,
            'kpoint_path': kpath_data,
        }
        wannier_lines = WannierGrammar().write(wannierdict).splitlines()

        # Add the lines to epwdict. 
        for line_num, line in enumerate(wannier_lines):
            epwdict['inputepw'][f'wdata({line_num+1})'] = f"'{line}'"

    @property
    def elph(self) -> str:
        epwdict: dict = {
            'inputepw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'dvscf_dir': "'./save'",
                
                'nk1': jmespath.search('elph.coarse_kgrid[0]', self.inputdict),
                'nk2': jmespath.search('elph.coarse_kgrid[1]', self.inputdict),
                'nk3': jmespath.search('elph.coarse_kgrid[2]', self.inputdict),
                'nq1': jmespath.search('elph.coarse_qgrid[0]', self.inputdict),
                'nq2': jmespath.search('elph.coarse_qgrid[1]', self.inputdict),
                'nq3': jmespath.search('elph.coarse_qgrid[2]', self.inputdict),
                'nkf1': jmespath.search('elph.fine_kgrid[0]', self.inputdict),
                'nkf2': jmespath.search('elph.fine_kgrid[1]', self.inputdict),
                'nkf3': jmespath.search('elph.fine_kgrid[2]', self.inputdict),
                'nqf1': jmespath.search('elph.fine_qgrid[0]', self.inputdict),
                'nqf2': jmespath.search('elph.fine_qgrid[1]', self.inputdict),
                'nqf3': jmespath.search('elph.fine_qgrid[2]', self.inputdict),
                'nbndsub': jmespath.search('elph.val_bands', self.inputdict) + jmespath.search('elph.cond_bands', self.inputdict),

                
                'elph': '.true.',
                'epbwrite': '.true.',
                'epbread': '.false.',
                'epwwrite': '.true.',
                'epwread': '.false.',
                'lpolar': '.true.',

                # # Reuse wannier90 calculations. 
                # 'wannierize': '.false.',
                # 'filukk': "'./struct.ukk'",

                # Do custom wannier calculations. 
                'wannierize': '.true.',
                'wannier_plot': '.true.',
                'wannier_plot_supercell': f'{jmespath.search("elph.coarse_kgrid[0]", self.inputdict)*2} {jmespath.search("elph.coarse_kgrid[1]", self.inputdict)*2} {jmespath.search("elph.coarse_kgrid[2]", self.inputdict)*2}',
                'auto_projections': '.true.',
                'scdm_proj': '.true.',
            }
        }

        # Add bands_skipped.
        a = self.bands_skipped_string
        if self.bands_skipped_string is not None and self.bands_skipped_string!='' : 
            epwdict['inputepw']['bands_skipped'] = self.bands_skipped_string

        # Add extra wannier data. 
        self.add_wannier_data(epwdict)

        # Update if needed. 
        update_dict(epwdict, jmespath.search('elph.args', self.inputdict))

        # If projections provided, delete auto_projections and scdm_proj.
        if 'proj(1)' in epwdict['inputepw'].keys() and epwdict['inputepw']['wannierize'] == '.true.': 
            del epwdict['inputepw']['auto_projections']
            del epwdict['inputepw']['scdm_proj']

        return NamelistGrammar().write(epwdict)

    @property
    def job_elph(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'elph.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

rm -rf ./tmp
cp -r ../{jmespath.search('elph.nscf_link', self.inputdict)}/tmp ./tmp
rm -rf ./save
cp -r ../dfpt/save ./save
{scheduler.get_exec_prefix()}epw.x {scheduler.get_exec_infix()} < elph.in  &> elph.in.out 
cp ./tmp/struct.xml ./save/wfn.xml
cp ./tmp/*epb* ./save/

wannierqe_pp="
from fpflow.analysis.wannierqe import WannierQeAnalysis
wan = WannierQeAnalysis()
wan.read_all()
wan.write()
"

python -c "$wannierqe_pp" &> wannierqe_pp.out 
'''
        return file_string


    @property
    def file_contents(self) -> dict:
        return {
            './elph_grid/elph.in': self.elph,
            './elph_grid/job_elph.sh': self.job_elph,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './elph_grid/job_elph.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './elph_grid',
        ]
    
    def plot(self, **kwargs):
            EpwElphGridPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion