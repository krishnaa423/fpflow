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
from fpflow.plots.phbands import PhbandsPlot
from fpflow.structure.kpath import Kpath
from fpflow.inputs.grammars.wannier import WannierGrammar
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.structure.kpts import Kpts
import xml.etree.ElementTree as ET
from ase.units import Hartree
import numpy as np 
from fpflow.io.read_write import str_2_f
from fpflow.plots.wannier import QeWannierPlot

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeWannierStep(Step):
    def add_exclude_bands(self, wannierdict: dict, max_val_bands: int, val_bands: int):
        if val_bands < max_val_bands:
            if val_bands == max_val_bands -1:
                wannierdict['exclude_bands'] = ' 1 '
            else:
                wannierdict['exclude_bands'] = f' 1:{max_val_bands - val_bands} '

    @property
    def wan(self) -> str:
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = int(qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))
        cond_bands: int = jmespath.search('wannier.cond_bands', self.inputdict)
        val_bands: int = jmespath.search('wannier.val_bands', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('wannier.kgrid[0]', self.inputdict),
                jmespath.search('wannier.kgrid[1]', self.inputdict),
                jmespath.search('wannier.kgrid[2]', self.inputdict),
            ],
            is_reduced=False,
        ).wannier_kpts

        # Kpath.
        kpath_data: list = Kpath.from_yamlfile().wannierpath_list

        self.wannierdict: dict = {
            'mp_grid': jmespath.search('wannier.kgrid', self.inputdict),
            'num_bands': cond_bands + val_bands,
            'num_wann': cond_bands + val_bands,
            'auto_projections': True,
            'wannier_plot': True,
            'bands_plot': True,
            'bands_num_points': jmespath.search('kpath.npoints_segment', self.inputdict),
            'write_hr': True,
            'use_ws_distance': True,
            'write_xyz': True,
            'write_u_matrices': True,
            'write_tb': True,
            'write_bvec': True,
            'unit_cell_cart': qestruct.cell,
            'atoms_frac': qestruct.scaled_atomic_positions,
            'kpoints': kpts,
            'kpoint_path': kpath_data,
        }

        if jmespath.search('scf.is_spinorbit', self.inputdict):
            self.wannierdict['spinors'] = True

        # Add exclude bands if needed.
        self.add_exclude_bands(self.wannierdict, max_val_bands, val_bands)

        # Update if needed. 
        update_dict(self.wannierdict, jmespath.search('wannier.args', self.inputdict))

        # If projections provided, delete auto_projections.
        if 'projections' in self.wannierdict.keys(): del self.wannierdict['auto_projections']

        return WannierGrammar().write(self.wannierdict)
    
    @property
    def pw2wan(self) -> str:
        pw2wandict: dict = {
            'inputpp': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'seedname': "'struct'",
                'write_amn': '.true.',
                'write_mmn': '.true.',
                'write_unk': '.true.',
                'scdm_proj': '.true.',
            }
        }

        # Update if needed. 
        update_dict(pw2wandict, jmespath.search('wannier.pw2wan_args', self.inputdict))

        # If projections provided, delete scdm_proj.
        if 'projections' in self.wannierdict: del pw2wandict['inputpp']['scdm_proj']

        file_string: str =  NamelistGrammar().write(pw2wandict)

        return file_string
    
    @property
    def create_ukk(self) -> str:
        filestring = '''
using WannierIO

# Read wannier90 files
chk = WannierIO.read_chk("./struct.chk")

# Read qe xml.
qe_xml = WannierIO.read_qe_xml("./tmp/struct.xml")
alat = qe_xml.alat

# Create ukk file.
ukk = WannierIO.Ukk(chk, alat)
WannierIO.write_epw_ukk("struct.ukk", ukk)

println("Created struct.ukk file for EPW calculations.")
'''

        return filestring

    @property
    def job_wanpp(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wannier.job_wanpp_info')
        
        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf ../{jmespath.search('wannier.nscf_link', self.inputdict)}/tmp ./tmp
{scheduler.get_exec_prefix()}wannier90.x {scheduler.get_exec_infix()} -pp struct &> struct.win.pp.out
'''
        return file_string
    
    @property
    def job_pw2wan(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wannier.job_pw2wan_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw2wannier90.x {scheduler.get_exec_infix()} < pw2wan.in &> pw2wan.in.out
'''
        return file_string
    
    @property
    def job_wan(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'wannier.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}wannier90.x {scheduler.get_exec_infix()} struct &> struct.win.out

# Create .ukk file for epw.
# julia create_ukk.jl &> create_ukk.jl.out

wannierqe_pp="
from fpflow.analysis.wannierqe import WannierQeAnalysis
wan = WannierQeAnalysis()
wan.read_all()
wan.write()
"

tbmodels_pp='
import tbmodels

model = tbmodels.Model.from_wannier_files(
    hr_file="./struct_hr.dat",
    wsvec_file="./struct_wsvec.dat",
    xyz_file="./struct_centres.xyz",
    win_file="./struct.win",
)
model.to_hdf5_file("./wannier_tbmodel.hdf5")
'

python -c "$wannierqe_pp" &> wannierqe_pp.out 
python -c "$tbmodels_pp" &> wannierqe_pp_tbmodels.out 
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        file_list = {
            './wannier/struct.win': self.wan,
            './wannier/pw2wan.in': self.pw2wan,
            # './wannier/create_ukk.jl': self.create_ukk,
            './wannier/job_wanpp.sh': self.job_wanpp,
            './wannier/job_pw2wan.sh': self.job_pw2wan,
            './wannier/job_wan.sh': self.job_wan,
        }

        return file_list
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './wannier/job_wanpp.sh',
            './wannier/job_pw2wan.sh',
            './wannier/job_wan.sh',
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './wannier',
        ]

    def plot(self, **kwargs):
            QeWannierPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion