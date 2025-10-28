#region modules
from typing import List
import jmespath
from webcolors import names 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from benedict import benedict
from fpflow.inputs.grammars.namelist import NamelistGrammar
from fpflow.inputs.grammars.qe import QeGrammar
from fpflow.inputs.grammars.bgw import BgwGrammar
from fpflow.io.update import update_dict
from fpflow.schedulers.scheduler import Scheduler
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.structure.kpts import Kpts

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeNscfStep(Step):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_wfn(self, options: dict) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = int(qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))
        cond_bands: int = jmespath.search('cond_bands', options)

        # Kpts.
        qshift = jmespath.search('qshift', options)
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('kgrid[0]', options),
                jmespath.search('kgrid[1]', options),
                jmespath.search('kgrid[2]', options),
            ],
            qshift=qshift,
            is_reduced=jmespath.search('sym', options),
        )

        wfndict: dict = {
            'control': {
                'outdir': './tmp',
                'prefix': 'struct',
                'pseudo_dir': './pseudos',
                'calculation': 'bands',
                'tprnfor': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                'nbnd': int(cond_bands + max_val_bands),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict)
            },
            'electrons': {},
            'ions': {},
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'crystal',
                'nkpt': kpts.nkpt,
                'data': kpts.wfn_kpts if qshift==[0.0, 0.0, 0.0] else kpts.wfnq_kpts,
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            wfndict['system']['noncolin'] = True
            wfndict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(wfndict, jmespath.search('args', options))

        return QeGrammar().write(wfndict)

    def get_pw2bgw(self, options: dict) -> str:
        # Qestruct.
        max_val_bands: int = int(QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))
        cond_bands: int = jmespath.search('cond_bands', options)

        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./tmp'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'wfn'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('kgrid[0]', options),
                'wfng_nk2': jmespath.search('kgrid[1]', options),
                'wfng_nk3': jmespath.search('kgrid[2]', options),
                'wfng_dk1': jmespath.search('qshift[0]', options) * jmespath.search('kgrid[0]', options),
                'wfng_dk2': jmespath.search('qshift[1]', options) * jmespath.search('kgrid[1]', options),
                'wfng_dk3': jmespath.search('qshift[2]', options) * jmespath.search('kgrid[2]', options),
            }
        }

        # Extract rho, vxc, vsc, vkb if asked.
        use_rho_vxc_vsc_vkb: bool = jmespath.search('extract_rho_vxc_vsc_vkb', options)
        if use_rho_vxc_vsc_vkb:
            pw2bgwdict['input_pw2bgw'].update({
                'rhog_flag': '.true.',
                'rhog_file': "'rho'",
                'vxc_flag': '.true.',
                'vxc_file': "'vxc.dat'",
                'vxc_diag_nmin': 1,
                'vxc_diag_nmax': max_val_bands + cond_bands,
                'vxc_offdiag_nmin': 0,
                'vxc_offdiag_nmax': 0,
                'vscg_flag': '.true.',
                'vscg_file': "'vsc'",
                'vkbg_flag': '.true.',
                'vkbg_file': "'vkb'",
            })

        # Update from yaml if needed. 
        update_dict(pw2bgwdict, jmespath.search('pw2bgw_args', options))

        return NamelistGrammar().write(pw2bgwdict)

    def get_parabands(self, options: dict) -> str:
        max_val_bands: int = QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )
        parabands_cond_bands: int = jmespath.search('parabands.cond_bands', options)

        parabandsdict: dict = {
            'input_wfn_file': 'wfn',
            'output_wfn_file': 'wfn_parabands.h5',
            'vsc_file': 'vsc',
            'vkb_file': 'vkb',
            'number_bands': parabands_cond_bands + max_val_bands,
            'wfn_io_mpiio_mode': 1,
        }

        # Add pseudobands option if needed.
        if jmespath.search('parabands.use_pseudobands', options):
            parabandsdict['use_pseudobands'] = ''

        # Update if needed. 
        update_dict(parabandsdict, jmespath.search('parabands.parabands_args', options))

        return BgwGrammar().write(parabandsdict)

    def get_job_wfn(self, options: dict) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(options, 'job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

rm -rf ./tmp
cp -r ../scf/tmp ./tmp
{scheduler.get_exec_prefix()}pw.x {scheduler.get_exec_infix()} < wfn.in &> wfn.in.out 

cp ./tmp/struct.xml ./wfn.xml
'''
        return file_string

    def get_job_pw2bgw(self, options: dict) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(options, 'job_pw2bgw_info')

        copy_extra_string: str = '''# Copy extra files if needed.
cp ./tmp/rho ./
cp ./tmp/vxc.dat ./
cp ./tmp/vsc ./
cp ./tmp/vkb ./'''

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}pw2bgw.x -pd .true. < pw2bgw.in &> pw2bgw.in.out
cp ./tmp/wfn ./
wfn2hdf.x BIN wfn wfn.h5

{copy_extra_string if jmespath.search('extract_rho_vxc_vsc_vkb', options) else ""}
'''
        return file_string

    def get_job_parabands(self, options: dict) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(options, 'parabands.job_parabands_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

{scheduler.get_exec_prefix()}parabands.cplx.x &> parabands.inp.out
'''
        return file_string

    @property
    def file_contents(self) -> dict:
        contents_dict: dict = {}

        for wfn_options in jmespath.search('nscf.list[*]', self.inputdict):
            name = jmespath.search('name', wfn_options)

            # Check if name is in enabled dirs.
            enabled_dirs: List[str] = jmespath.search('nscf.enabled_dirs', self.inputdict) or []
            if name in enabled_dirs:
                contents_dict.update({
                    f'./{name}/wfn.in': self.get_wfn(wfn_options),
                    f'./{name}/job_wfn.sh': self.get_job_wfn(wfn_options),
                    f'./{name}/pw2bgw.in': self.get_pw2bgw(wfn_options),
                    f'./{name}/job_pw2bgw.sh': self.get_job_pw2bgw(wfn_options),
                })

                # Add parabands if asked.
                use_parabands: bool = jmespath.search('parabands.enabled', wfn_options)
                if use_parabands:
                    contents_dict.update({
                        f'./{name}/parabands.inp': self.get_parabands(wfn_options),
                        f'./{name}/job_parabands.sh': self.get_job_parabands(wfn_options),
                    })

        return contents_dict

    @property
    def job_scripts(self) -> List[str]:
        scripts_list: list[str] = []

        for wfn_options in jmespath.search('nscf.list[*]', self.inputdict):
            name = jmespath.search('name', wfn_options)

            # Check if name is in enabled dirs.
            enabled_dirs: List[str] = jmespath.search('nscf.enabled_dirs', self.inputdict) or []
            if name in enabled_dirs:
                scripts_list.extend([
                    f'./{name}/job_wfn.sh',
                    f'./{name}/job_pw2bgw.sh',
                ])

                # Add parabands if asked.
                use_parabands: bool = jmespath.search('parabands.enabled', wfn_options)
                if use_parabands:
                    scripts_list.append(f'./{name}/job_parabands.sh')

        return scripts_list

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        inodes_list: List[str] = [f'./{dirname}' for dirname in jmespath.search('nscf.list[*].name', self.inputdict) if dirname in jmespath.search('nscf.enabled_dirs', self.inputdict) or []]

        return inodes_list

#endregion