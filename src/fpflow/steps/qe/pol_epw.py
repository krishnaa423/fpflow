#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.structure.qe.qe_struct import QeStruct
import jmespath
from fpflow.io.update import update_dict
from fpflow.inputs.grammars.namelist import NamelistGrammar
from fpflow.schedulers.scheduler import Scheduler
from fpflow.plots.pol_epw import EpwPolPlot

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class EpwPolStep(Step):
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

    def get_added_dict_str(self, added_dict: dict) -> dict:
        # wfnlink. 
        wfn_link: str = jmespath.search('elph.nscf_link', self.inputdict)
        wfn_names: list[str] = jmespath.search(f'nscf.list[*].name', self.inputdict)
        wfn_index: int = wfn_names.index(wfn_link)
        self.wfn_options: dict = jmespath.search(f'nscf.list[{wfn_index}]', self.inputdict)

        elphdict: dict = {
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
                'nbndsub': jmespath.search('elph.val_bands', self.inputdict) + jmespath.search('elph.cond_bands', self.inputdict),

                'lpolar': '.true.',
            }
        }

        # Add bands_skipped.
        a = self.bands_skipped_string
        if self.bands_skipped_string is not None and self.bands_skipped_string!='' : 
            elphdict['inputepw']['bands_skipped'] = self.bands_skipped_string

        # Update if needed. 
        update_dict(elphdict, added_dict)
        update_dict(elphdict, jmespath.search('pol_epw.args', self.inputdict))

        return NamelistGrammar().write(elphdict)

    @property
    def elph_in(self) -> str:
        added_dict: dict = {
            'inputepw': {
                'elph': '.true.',
                'wannierize': '.true.',
                'wannier_plot':'.true.',
                'wannier_plot_supercell': f"{jmespath.search('xctph.qgrid[0]', self.inputdict)*2} {jmespath.search('xctph.qgrid[1]', self.inputdict)*2} {jmespath.search('xctph.qgrid[2]', self.inputdict)*2}",
                'auto_projections': '.true.',
                'scdm_proj': '.true.',
                'iprint': 2,
            }
        }

        return self.get_added_dict_str(added_dict)

    @property
    def job_elph(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'pol.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

rm -rf ./tmp
ln -sf ../{self.wfn_options["name"]}/tmp ./tmp
ln -sf ../dfpt/save ./save
{scheduler.get_exec_prefix()}epw.x -npool {scheduler.nk} < elph.in  &> elph.in.out
'''
        return file_string

    @property
    def elec_in(self) -> str:
        added_dict: dict = {
            'inputepw': {
                'elph': '.true.',
                'epwread': '.true.',
                'wannierize': '.false.',
                
                'plrn': '.true.',
                'restart_plrn': '.false.',
                'type_plrn': -1,    # Electron
                'init_plrn': 1,  # Gaussian initial guess.
                'init_sigma_plrn': 1.0,  # Broadening for the gaussian,
                'niter_plrn': jmespath.search('xctpol.max_steps', self.inputdict),
                'conv_thr_plrn': jmespath.search('xctpol.max_error', self.inputdict),
                'ethrdg_plrn': '1e-5',
            }
        }

        return self.get_added_dict_str(added_dict)

    @property
    def job_elec(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'pol.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}epw.x -npool {scheduler.nk} < elec.in  &> elec.in.out
'''
        return file_string

    @property
    def elec_pp_in(self) -> str:
        added_dict: dict = {
            'inputepw': {
                'elph': '.true.',
                'epwread': '.true.',
                'wannierize': '.false.',
                
                'plrn': '.true.',
                'restart_plrn': '.true.',
                'type_plrn': -1,    # Electron
                'cal_psir_plrn' : '.true.',
                'interp_Ank_plrn': '.true.',
                'interp_Bqu_plrn': '.true.',
            }
        }

        return self.get_added_dict_str(added_dict)

    @property
    def job_elec_pp(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'pol.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}epw.x -npool {scheduler.nk} < elec_pp.in  &> elec_pp.in.out

# Move files.
mv Amp.plrn elec_Amp.plrn
mv Ank.band.plrn elec_Ank.band.plrn
mv Ank.plrn elec_Ank.plrn
mv Bmat.plrn elec_Bmat.plrn
mv Bmat.band.plrn elec_Bmat.band.plrn
mv dos.plrn elec_dos.plrn
mv dtau.plrn elec_dtau.plrn
mv dtau.plrn.xsf elec_dtau.plrn.xsf
mv psir_plrn.xsf elec_psir_plrn.xsf
'''
        return file_string

    @property
    def hole_in(self) -> str:
        added_dict: dict = {
            'inputepw': {
                'elph': '.true.',
                'epwread': '.true.',
                'wannierize': '.false.',
                
                'plrn': '.true.',
                'restart_plrn': '.false.',
                'type_plrn': 1,    # Hole
                'init_plrn': 1,  # Gaussian initial guess.
                'init_sigma_plrn': 1.0,  # Broadening for the gaussian,
                'niter_plrn': jmespath.search('xctpol.max_steps', self.inputdict),
                'conv_thr_plrn': jmespath.search('xctpol.max_error', self.inputdict),
                'ethrdg_plrn': '1e-5',
            }
        }

        return self.get_added_dict_str(added_dict)

    @property
    def job_hole(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'pol.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}epw.x -npool {scheduler.nk} < hole.in  &> hole.in.out
'''
        return file_string

    @property
    def hole_pp_in(self) -> str:
        added_dict: dict = {
            'inputepw': {
                'elph': '.true.',
                'epwread': '.true.',
                'wannierize': '.false.',
                
                'plrn': '.true.',
                'restart_plrn': '.true.',
                'type_plrn': 1,    # Hole
                'cal_psir_plrn' : '.true.',
                'interp_Ank_plrn': '.true.',
                'interp_Bqu_plrn': '.true.',
            }
        }

        return self.get_added_dict_str(added_dict)

    @property
    def job_hole_pp(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'pol.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}epw.x -npool {scheduler.nk} < hole_pp.in  &> hole_pp.in.out

# Move files.
mv Amp.plrn hole_Amp.plrn
mv Ank.band.plrn hole_Ank.band.plrn
mv Ank.plrn hole_Ank.plrn
mv Bmat.plrn hole_Bmat.plrn
mv Bmat.band.plrn hole_Bmat.band.plrn
mv dos.plrn hole_dos.plrn
mv dtau.plrn hole_dtau.plrn
mv dtau.plrn.xsf hole_dtau.plrn.xsf
mv psir_plrn.xsf hole_psir_plrn.xsf
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        return {
            './pol_epw/elph.in': self.elph_in,
            './pol_epw/job_elph.sh': self.job_elph,
            './pol_epw/elec.in': self.elec_in,
            './pol_epw/job_elec.sh': self.job_elec,
            './pol_epw/elec_pp.in': self.elec_pp_in,
            './pol_epw/job_elec_pp.sh': self.job_elec_pp,
            './pol_epw/hole.in': self.hole_in,
            './pol_epw/job_hole.sh': self.job_hole,
            './pol_epw/hole_pp.in': self.hole_pp_in,
            './pol_epw/job_hole_pp.sh': self.job_hole_pp,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './pol_epw/job_elph.sh',
            './pol_epw/job_elec.sh',
            './pol_epw/job_elec_pp.sh',
            './pol_epw/job_hole.sh',
            './pol_epw/job_hole_pp.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './pol_epw',
        ]

    def plot(self, **kwargs):
        EpwPolPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion