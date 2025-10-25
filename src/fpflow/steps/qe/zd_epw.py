#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from importlib.util import find_spec 
from fpflow.inputs.grammars.namelist import NamelistGrammar
from fpflow.inputs.grammars.qe import QeGrammar
from fpflow.inputs.grammars.bgw import BgwGrammar
import jmespath
from fpflow.io.update import update_dict
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.schedulers.scheduler import Scheduler
from fpflow.structure.kpts import Kpts
from fpflow.plots.zd_epw import EpwZdPlot

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class EpwZdStep(Step):
    @property
    def setup(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'zd.job_info')
        scheduler_epw: Scheduler = Scheduler.from_jmespath(self.inputdict, 'zd.job_epw_info')

        pool:int = scheduler_epw.nk if scheduler_epw.nk is not None else 1

        file_string = f'''#!/bin/bash

MPIRUN="{scheduler.mpi_exec} -n {scheduler.ntasks} "
MPIRUN_EPW="{scheduler_epw.mpi_exec} -n {pool} "
NQPT={jmespath.search('common.fine_nk', self.inputdict)}
NPOOL={pool}
'''
        
        return file_string

    @property
    def espresso_scf(self) -> str:
        qestruct = QeStruct.from_inputdict(self.inputdict)

        qedict: dict = {
            'control': {
                'outdir': './',
                'prefix': 'struct',
                'pseudo_dir': '../../pseudos/',
                'calculation': 'scf',
                'tprnfor': True,
                'tstress': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict),
                'degauss': 1e-8,
            },
            'electrons': {
                'conv_thr': 1e-12,
            },
            'ions': {},
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'automatic',
                'data': [
                    jmespath.search('scf.kgrid[0]', self.inputdict),
                    jmespath.search('scf.kgrid[1]', self.inputdict),
                    jmespath.search('scf.kgrid[2]', self.inputdict),
                    0,
                    0,
                    0,
                ],
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            qedict['system']['noncolin'] = True
            qedict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(qedict, jmespath.search('scf.args', self.inputdict))

        return QeGrammar().write(qedict)

    @property
    def espresso_wfn(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = int(qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))
        cond_bands: int = jmespath.search('wfn.cond_bands', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('wfn.kgrid[0]', self.inputdict),
                jmespath.search('wfn.kgrid[1]', self.inputdict),
                jmespath.search('wfn.kgrid[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfn.sym', self.inputdict),
        )

        wfndict: dict = {
            'control': {
                'outdir': './',
                'prefix': 'struct',
                'pseudo_dir': '../../pseudos/',
                'calculation': 'bands',
                'tprnfor': True,
                'tstress': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                'nbnd': int(cond_bands + max_val_bands),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict),
                'degauss': 1e-8,
            },
            'electrons': {
                'conv_thr': 1e-12,
                'diagonalization': 'cg',
                'diago_full_acc': True,
            },
            'ions': {},
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'crystal',
                'nkpt': kpts.nkpt,
                'data': kpts.wfn_kpts,
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            wfndict['system']['noncolin'] = True
            wfndict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(wfndict, jmespath.search('wfn.args', self.inputdict))

        return QeGrammar().write(wfndict)

    @property
    def espresso_wfn_pp(self) -> str:
        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'wfn.cplx'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('wfn.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('wfn.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('wfn.kgrid[2]', self.inputdict),
                'wfng_dk1': 0.0,
                'wfng_dk2': 0.0,
                'wfng_dk3': 0.0,
                'vscg_flag': '.true.',
                'vscg_file': "'VSC'",
                'vkbg_flag': '.true.',
                'vkbg_file': "'VKB'",
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('wfn.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def espresso_wfn_parabands(self) -> str:
        max_val_bands: int = QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )
        parabands_cond_bands: int = jmespath.search('wfn.parabands_cond_bands', self.inputdict)

        parabandsdict: dict = {
            'input_wfn_file': 'wfn.cplx',
            'output_wfn_file': 'wfn.cplx.h5',
            'vsc_file': 'VSC',
            'vkb_file': 'VKB',
            'number_bands': parabands_cond_bands + max_val_bands,
            'wfn_io_mpiio_mode': 1,
        }

        # Add pseudobands option if needed.
        if jmespath.search('wfn.is_pseudobands', self.inputdict):
            parabandsdict['use_pseudobands'] = ''

        # Update if needed. 
        update_dict(parabandsdict, jmespath.search('wfn.parabands_args', self.inputdict))

        return BgwGrammar().write(parabandsdict)

    @property
    def espresso_wfnq(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('wfnq.kgrid[0]', self.inputdict),
                jmespath.search('wfnq.kgrid[1]', self.inputdict),
                jmespath.search('wfnq.kgrid[2]', self.inputdict),
            ],
            qshift=[
                jmespath.search('wfnq.qshift[0]', self.inputdict),
                jmespath.search('wfnq.qshift[1]', self.inputdict),
                jmespath.search('wfnq.qshift[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfnq.sym', self.inputdict),
        )

        wfndict: dict = {
            'control': {
                'outdir': './',
                'prefix': 'struct',
                'pseudo_dir': '../../pseudos/',
                'calculation': 'bands',
                'tprnfor': True,
                'tstress': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                # 'nbnd': int(cond_bands + max_val_bands),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict),
                'degauss': 1e-8,
            },
            'electrons': {
                'conv_thr': 1e-12,
                'diago_full_acc': True,
            },
            'ions': {},
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'crystal',
                'nkpt': kpts.nkpt,
                'data': kpts.wfnq_kpts,
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            wfndict['system']['noncolin'] = True
            wfndict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(wfndict, jmespath.search('wfnq.args', self.inputdict))

        return QeGrammar().write(wfndict)

    @property
    def espresso_wfnq_pp(self) -> str:
        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'wfn.cplx'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('wfnq.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('wfnq.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('wfnq.kgrid[2]', self.inputdict),
                'wfng_dk1': jmespath.search('wfnq.qshift[0]', self.inputdict) * jmespath.search('wfnq.kgrid[0]', self.inputdict),
                'wfng_dk2': jmespath.search('wfnq.qshift[1]', self.inputdict) * jmespath.search('wfnq.kgrid[1]', self.inputdict),
                'wfng_dk3': jmespath.search('wfnq.qshift[2]', self.inputdict) * jmespath.search('wfnq.kgrid[2]', self.inputdict),
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('wfnq.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def espresso_wfnco(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = int(qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))
        cond_bands: int = jmespath.search('wfnfi.cond_bands', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('wfnfi.kgrid[0]', self.inputdict),
                jmespath.search('wfnfi.kgrid[1]', self.inputdict),
                jmespath.search('wfnfi.kgrid[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfnfi.sym', self.inputdict),
        )

        wfndict: dict = {
            'control': {
                'outdir': './',
                'prefix': 'struct',
                'pseudo_dir': '../../pseudos/',
                'calculation': 'bands',
                'tprnfor': True,
                'tstress': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                'nbnd': int(cond_bands + max_val_bands),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict),
                'degauss': 1e-8,
            },
            'electrons': {
                'conv_thr': 1e-12,
                'diagonalization': 'cg',
                'diago_full_acc': True,
            },
            'ions': {},
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'crystal',
                'nkpt': kpts.nkpt,
                'data': kpts.wfn_kpts,
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            wfndict['system']['noncolin'] = True
            wfndict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(wfndict, jmespath.search('wfnfi.args', self.inputdict))

        return QeGrammar().write(wfndict)

    @property
    def espresso_wfnco_pp(self) -> str:
        # Qestruct.
        max_val_bands: int = int(QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))

        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'wfn.cplx'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('wfnfi.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('wfnfi.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('wfnfi.kgrid[2]', self.inputdict),
                'wfng_dk1': 0.0,
                'wfng_dk2': 0.0,
                'wfng_dk3': 0.0,
                'rhog_flag': '.true.',
                'rhog_file': "'rho.real'",
                'vxc_flag': '.true.',
                'vxc_file': "'vxc.dat'",
                'vxcg_flag': '.true.',
                'vxcg_file': "'vxc.real'",
                'vxc_diag_nmin': max_val_bands - jmespath.search('gw.sigma.val_bands', self.inputdict) + 1,
                'vxc_diag_nmax': max_val_bands + jmespath.search('gw.sigma.cond_bands', self.inputdict),
                'vscg_flag': '.true.',
                'vscg_file': "'VSC'",
                'vkbg_flag': '.true.',
                'vkbg_file': "'VKB'",
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('wfnfi.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def espresso_wfnco_parabands(self) -> str:
        max_val_bands: int = QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )
        parabands_cond_bands: int = jmespath.search('wfn.parabands_cond_bands', self.inputdict)

        parabandsdict: dict = {
            'input_wfn_file': 'wfn.cplx',
            'output_wfn_file': 'wfn.cplx.h5',
            'vsc_file': 'VSC',
            'vkb_file': 'VKB',
            'number_bands': parabands_cond_bands + max_val_bands,
            'wfn_io_mpiio_mode': 1,
        }

        # Add pseudobands option if needed.
        if jmespath.search('wfn.is_pseudobands', self.inputdict):
            parabandsdict['use_pseudobands'] = ''

        # Update if needed. 
        update_dict(parabandsdict, jmespath.search('wfn.parabands_args', self.inputdict))

        return BgwGrammar().write(parabandsdict)

    @property
    def espresso_wfnfi(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = int(qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))
        cond_bands: int = jmespath.search('wfnfi.cond_bands', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('wfnfi.kgrid[0]', self.inputdict),
                jmespath.search('wfnfi.kgrid[1]', self.inputdict),
                jmespath.search('wfnfi.kgrid[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfnfi.sym', self.inputdict),
        )

        wfndict: dict = {
            'control': {
                'outdir': './',
                'prefix': 'struct',
                'pseudo_dir': '../../pseudos/',
                'calculation': 'bands',
                'tprnfor': True,
                'tstress': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                'nbnd': int(cond_bands + max_val_bands),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict),
                'degauss': 1e-8,
            },
            'electrons': {
                'conv_thr': 1e-12,
                'diago_full_acc': True,
            },
            'ions': {},
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'crystal',
                'nkpt': kpts.nkpt,
                'data': kpts.wfn_kpts,
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            wfndict['system']['noncolin'] = True
            wfndict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(wfndict, jmespath.search('wfnfi.args', self.inputdict))

        return QeGrammar().write(wfndict)

    @property
    def espresso_wfnfi_pp(self) -> str:
         # Qestruct.
        max_val_bands: int = int(QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))

        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'wfn.cplx'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('wfnfi.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('wfnfi.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('wfnfi.kgrid[2]', self.inputdict),
                'wfng_dk1': 0.0,
                'wfng_dk2': 0.0,
                'wfng_dk3': 0.0,
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('wfnfi.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def espresso_wfnqfi(self) -> str:
        # Qestruct.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        max_val_bands: int = int(qestruct.max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))
        cond_bands: int = jmespath.search('wfnqfi.cond_bands', self.inputdict)

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('wfnqfi.kgrid[0]', self.inputdict),
                jmespath.search('wfnqfi.kgrid[1]', self.inputdict),
                jmespath.search('wfnqfi.kgrid[2]', self.inputdict),
            ],
            qshift=[
                jmespath.search('wfnq.qshift[0]', self.inputdict),
                jmespath.search('wfnq.qshift[1]', self.inputdict),
                jmespath.search('wfnq.qshift[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfnqfi.sym', self.inputdict),
        )

        wfndict: dict = {
            'control': {
                'outdir': './',
                'prefix': 'struct',
                'pseudo_dir': '../../pseudos/',
                'calculation': 'bands',
                'tprnfor': True,
                'tstress': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                # 'nbnd': int(cond_bands + max_val_bands),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict),
                'degauss': 1e-8,
            },
            'electrons': {
                'conv_thr': 1e-12,
                'diago_full_acc': True,
            },
            'ions': {},
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'crystal',
                'nkpt': kpts.nkpt,
                'data': kpts.wfnq_kpts,
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            wfndict['system']['noncolin'] = True
            wfndict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(wfndict, jmespath.search('wfnqfi.args', self.inputdict))

        return QeGrammar().write(wfndict)

    @property
    def espresso_wfnqfi_pp(self) -> str:
        pw2bgwdict: dict = {
            'input_pw2bgw': {
                'outdir': "'./'",
                'prefix': "'struct'",
                'real_or_complex': '2',
                'wfng_flag': '.true.',
                'wfng_file': "'wfn.cplx'",
                'wfng_kgrid': '.true.',
                'wfng_nk1': jmespath.search('wfnqfi.kgrid[0]', self.inputdict),
                'wfng_nk2': jmespath.search('wfnqfi.kgrid[1]', self.inputdict),
                'wfng_nk3': jmespath.search('wfnqfi.kgrid[2]', self.inputdict),
                'wfng_dk1': jmespath.search('wfnq.qshift[0]', self.inputdict) * jmespath.search('wfnqfi.kgrid[0]', self.inputdict),
                'wfng_dk2': jmespath.search('wfnq.qshift[1]', self.inputdict) * jmespath.search('wfnqfi.kgrid[1]', self.inputdict),
                'wfng_dk3': jmespath.search('wfnq.qshift[2]', self.inputdict) * jmespath.search('wfnqfi.kgrid[2]', self.inputdict),
            }
        }

        # Update if needed. 
        update_dict(pw2bgwdict, jmespath.search('wfnqfi.pw2bgw_args', self.inputdict))

        return NamelistGrammar().write(pw2bgwdict)

    @property
    def epsilon(self) -> str:
        # Qestruct.
        max_val_bands: int = int(QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        ))

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('wfn.kgrid[0]', self.inputdict),
                jmespath.search('wfn.kgrid[1]', self.inputdict),
                jmespath.search('wfn.kgrid[2]', self.inputdict),
            ],
            qshift=[
                jmespath.search('wfnq.qshift[0]', self.inputdict),
                jmespath.search('wfnq.qshift[1]', self.inputdict),
                jmespath.search('wfnq.qshift[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfn.sym', self.inputdict),
        )

        epsilondict: dict = {
            'number_bands': jmespath.search('gw.epsilon.cond_bands', self.inputdict) + max_val_bands,
            'degeneracy_check_override': '',
            'epsilon_cutoff': jmespath.search('gw.epsilon.ecut', self.inputdict),
            # 'use_wfn_hdf5': '',
            'qpoints': kpts.epsilon_kpts,
            'use_wfn_hdf5': '',
        }

        # Update if needed. 
        update_dict(epsilondict, jmespath.search('gw.epsilon.args', self.inputdict))

        return BgwGrammar().write(epsilondict)

    @property
    def sigma(self) -> str:
        # Qestruct.
        max_val_bands: int = QeStruct.from_inputdict(self.inputdict).max_val(
            xc=jmespath.search('scf.xc', self.inputdict),
            is_soc=jmespath.search('scf.is_spinorbit', self.inputdict),
        )

        # Kpts.
        kpts: Kpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('wfnfi.kgrid[0]', self.inputdict),
                jmespath.search('wfnfi.kgrid[1]', self.inputdict),
                jmespath.search('wfnfi.kgrid[2]', self.inputdict),
            ],
            is_reduced=jmespath.search('wfnfi.sym', self.inputdict),
        )

        sigmadict: dict = {
            'number_bands': jmespath.search('gw.sigma.conv_cond_bands', self.inputdict) + max_val_bands,
            'band_index_min': max_val_bands - jmespath.search('gw.sigma.val_bands', self.inputdict) + 1,
            'band_index_max': max_val_bands + jmespath.search('gw.sigma.cond_bands', self.inputdict),
            'degeneracy_check_override': '',
            'no_symmetries_q_grid': '',
            'screened_coulomb_cutoff': jmespath.search('gw.sigma.ecut', self.inputdict),
            'kpoints': kpts.sigma_kpts,
            'use_wfn_hdf5': '',
        }

        # Update if needed. 
        update_dict(sigmadict, jmespath.search('gw.sigma.args', self.inputdict))

        return BgwGrammar().write(sigmadict)

    @property
    def kernel(self) -> str:
        qshift: list[int] = jmespath.search('bse.kernel.qshift[*]', self.inputdict)

        kerneldict: dict = {
            'exciton_Q_shift': f"2 {qshift[0]:15.10f} {qshift[1]:15.10f} {qshift[2]:15.10f}",
            'no_symmetries_coarse_grid': '',
            'number_val_bands': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_cond_bands': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'dont_check_norms': '',
            'use_wfn_hdf5': '',
        }

        # Update if needed. 
        update_dict(kerneldict, jmespath.search('bse.kernel.args', self.inputdict))

        return BgwGrammar().write(kerneldict)

    @property
    def absorption(self) -> str:
        qshift: list[int] = jmespath.search('bse.absorption.qshift[*]', self.inputdict)
        pol_dir: list[int] = jmespath.search('bse.absorption.pol_dir[*]', self.inputdict)

        absorptiondict: dict = {
            'exciton_Q_shift': f"2 {qshift[0]:15.10f} {qshift[1]:15.10f} {qshift[2]:15.10f}",
            'no_symmetries_coarse_grid': '',
            'no_symmetries_fine_grid': '',
            'no_symmetries_shifted_grid': '',
            'number_val_bands_coarse': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_val_bands_fine': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_cond_bands_coarse': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'number_cond_bands_fine': jmespath.search('bse.absorption.cond_bands', self.inputdict) - 1,
            'degeneracy_check_override': '',
            'diagonalization': '',
            # 'use_elpa': '',  
            'use_momentum': '',  
            # 'use_velocity': '',  
            'polarization': ' '.join(list(map(str, pol_dir))),
            'eqp_co_corrections': '',
            'dump_bse_hamiltonian': '',
            'energy_resolution': 0.1,
            'write_eigenvectors': jmespath.search('bse.absorption.nxct', self.inputdict),
            'dont_check_norms': '',
            'use_wfn_hdf5': '',
        }

        # Update if needed. 
        update_dict(absorptiondict, jmespath.search('bse.absorption.args', self.inputdict))

        return BgwGrammar().write(absorptiondict)

    @property
    def seed_kernel(self) -> str:
        qshift: list[int] = jmespath.search('bse.kernel.qshift[*]', self.inputdict)

        kerneldict: dict = {
            'no_symmetries_coarse_grid': '',
            'number_val_bands': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_cond_bands': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'dont_check_norms': '',
            'use_wfn_hdf5': '',
        }

        # Update if needed. 
        update_dict(kerneldict, jmespath.search('bse.kernel.args', self.inputdict))

        return BgwGrammar().write(kerneldict)

    @property
    def seed_absorption(self) -> str:
        qshift: list[int] = jmespath.search('bse.absorption.qshift[*]', self.inputdict)
        pol_dir: list[int] = jmespath.search('bse.absorption.pol_dir[*]', self.inputdict)

        absorptiondict: dict = {
            'no_symmetries_coarse_grid': '',
            'no_symmetries_fine_grid': '',
            'no_symmetries_shifted_grid': '',
            'number_val_bands_coarse': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_val_bands_fine': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_cond_bands_coarse': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'number_cond_bands_fine': jmespath.search('bse.absorption.cond_bands', self.inputdict) - 1,
            'degeneracy_check_override': '',
            'diagonalization': '',
            # 'use_elpa': '',  
            'use_momentum': '',  
            # 'use_velocity': '',  
            'polarization': ' '.join(list(map(str, pol_dir))),
            'eqp_co_corrections': '',
            'dump_bse_hamiltonian': '',
            'energy_resolution': 0.1,
            'write_eigenvectors': jmespath.search('bse.absorption.nxct', self.inputdict),
            'dont_check_norms': '',
            'use_wfn_hdf5': '',
        }

        # Update if needed. 
        update_dict(absorptiondict, jmespath.search('bse.absorption.args', self.inputdict))

        return BgwGrammar().write(absorptiondict)

    @property
    def epw_scf(self) -> str:
        qestruct = QeStruct.from_inputdict(self.inputdict)

        qedict: dict = {
            'control': {
                'outdir': './',
                'prefix': 'struct',
                'pseudo_dir': '../pseudos/',
                'calculation': 'scf',
                'tprnfor': True,
                'tstress': True,
            },
            'system': {
                'ibrav': 0,
                'ntyp': qestruct.ntyp(),
                'nat': qestruct.nat(),
                'ecutwfc': jmespath.search('scf.ecut', self.inputdict),
                'degauss': 1e-8,
            },
            'electrons': {
                'conv_thr': 1e-12,
            },
            'ions': {},
            'cell': {},
            'atomic_species': qestruct.atomic_species,
            'cell_parameters': qestruct.cell,
            'atomic_positions': qestruct.atomic_positions,
            'k_points': {
                'type': 'automatic',
                'data': [
                    jmespath.search('scf.kgrid[0]', self.inputdict),
                    jmespath.search('scf.kgrid[1]', self.inputdict),
                    jmespath.search('scf.kgrid[2]', self.inputdict),
                    0,
                    0,
                    0,
                ],
            }
        }
        if jmespath.search('scf.is_spinorbit', self.inputdict):
            qedict['system']['noncolin'] = True
            qedict['system']['lspinorb'] = True

        # Update if needed. 
        update_dict(qedict, jmespath.search('scf.args', self.inputdict))

        return QeGrammar().write(qedict)

    @property
    def epw_ph(self) -> str:
        dfptdict: dict = {
            'inputph': {
                'outdir': "'./'",
                'prefix': "'struct'",
                'epsil': '.true.',
                'zeu': '.true.',
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
            epw_val_bands = jmespath.search('epw.val_bands', self.inputdict)
            total_val_bands = max_val_bands
            epw_cond_bands = jmespath.search('epw.cond_bands', self.inputdict)
            wfn_cond = jmespath.search('wfn.cond_bands', self.inputdict)

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
    def epw_xctph(self) -> str:
        epwdict: dict = {
            'inputepw': {
                'prefix': "'struct'",
                'outdir': "'./'",

                'elph': '.true.',
                'epwwrite': '.true.',
                'epwread': '.false.',

                'lpolar': '.true.',
                'nbndsub': jmespath.search('xctph.val_bands', self.inputdict) + jmespath.search('xctph.cond_bands', self.inputdict),
                'dvscf_dir': "'./save'",

                'wannierize': '.true.',
                'num_iter': 500,
                'iprint': 2,

                'exciton': '.true.',
                'explrn': '.false.',

                'negnv_explrn': jmespath.search('xctph.nxct', self.inputdict),
                'nbndv_explrn': jmespath.search('xctph.val_bands', self.inputdict),
                'nbndc_explrn': jmespath.search('xctph.cond_bands', self.inputdict),

                'wannier_plot':'.true.',
                'wannier_plot_supercell': f"{jmespath.search('xctph.qgrid[0]', self.inputdict)*2} {jmespath.search('xctph.qgrid[1]', self.inputdict)*2} {jmespath.search('xctph.qgrid[2]', self.inputdict)*2}",
                'wannier_plot_radius': '3.5d0',
                
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

                'band_plot': '.false',
            }
        }

        # Add bands_skipped.
        a = self.bands_skipped_string
        if self.bands_skipped_string is not None and self.bands_skipped_string!='' : 
            epwdict['inputepw']['bands_skipped'] = self.bands_skipped_string

        # Update if needed. 
        update_dict(epwdict, jmespath.search('zd.args_xctph', self.inputdict))

        return NamelistGrammar().write(epwdict)

    @property
    def epw_ste(self) -> str:
        epwdict: dict = {
            'inputepw': {
                'prefix': "'struct'",
                'outdir': "'./'",
                
                'elph': '.true.',
                'epwwrite': '.false.',
                'epwread': '.true.',

                'lpolar': '.true.',
                'nbndsub': jmespath.search('ste.val_bands', self.inputdict) + jmespath.search('ste.cond_bands', self.inputdict),
                'dvscf_dir': "'./save'",

                'wannierize': '.false.',

                'exciton': '.true.',
                'explrn': '.true.',

                'negnv_explrn': jmespath.search('ste.nxct', self.inputdict),
                'nbndv_explrn': jmespath.search('ste.val_bands', self.inputdict),
                'nbndc_explrn': jmespath.search('ste.cond_bands', self.inputdict),
                
                'ethrdg_plrn': jmespath.search('ste.max_error', self.inputdict),
                'init_plrn': 5,
                'niter_plrn': jmespath.search('ste.max_steps', self.inputdict),

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
                
            }
        }

        # Add bands_skipped.
        a = self.bands_skipped_string
        if self.bands_skipped_string is not None and self.bands_skipped_string!='' : 
            epwdict['inputepw']['bands_skipped'] = self.bands_skipped_string

        # Update if needed. 
        update_dict(epwdict, jmespath.search('zd.args_ste', self.inputdict))

        return NamelistGrammar().write(epwdict)

    @property
    def epw_ste_eh_densities(self) -> str:
        epwdict: dict = {
            'inputepw': {
                'prefix': "'struct'",
                'outdir': "'./'",
                
                'elph': '.true.',
                'epwwrite': '.false.',
                'epwread': '.true.',

                'lpolar': '.true.',
                'nbndsub': jmespath.search('ste.val_bands', self.inputdict) + jmespath.search('ste.cond_bands', self.inputdict),
                'dvscf_dir': "'./save'",

                'exciton': '.true.',
                'explrn': '.true.',
                'plot_explrn_e': '.true.',
                'plot_explrn_h': '.true.',

                'negnv_explrn': jmespath.search('ste.nxct', self.inputdict),
                'nbndv_explrn': jmespath.search('ste.val_bands', self.inputdict),
                'nbndc_explrn': jmespath.search('ste.cond_bands', self.inputdict),
                
                'step_wf_grid_plrn': 2,

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
                
            }
        }

        # Add bands_skipped.
        a = self.bands_skipped_string
        if self.bands_skipped_string is not None and self.bands_skipped_string!='' : 
            epwdict['inputepw']['bands_skipped'] = self.bands_skipped_string

        # Update if needed. 
        update_dict(epwdict, jmespath.search('zd.args_ste', self.inputdict))

        return NamelistGrammar().write(epwdict)

    @property
    def job_zd(self) -> str:
        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'zd.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

# Change to zd and run. Then change back.
cd ./zd
./run_all.sh
cd ../


'''
        
        return file_string

    @property
    def file_contents(self) -> dict:

        # Copy directory. 
        os.system('mkdir -p ./zd')
        zd_dir = os.path.join(*[
            os.path.dirname(find_spec('fpflow').origin),
            'data',
            'zd'
        ])
        os.system(f'cp -r {zd_dir}/* ./zd/')

        # Copy pseudos.
        qestruct = QeStruct.from_inputdict(self.inputdict)
        xc = jmespath.search('scf.xc', self.inputdict)
        is_soc = jmespath.search('scf.is_spinorbit', self.inputdict)

        paths, filenames = qestruct.get_pseudos_list(xc=xc, is_soc=is_soc)

        os.system('mkdir -p ./zd/pseudos/')
        for path, filename in zip(paths, filenames):
            os.system(f'cp {path} ./zd/pseudos/{filename}')

        return {
            './zd/setup.sh': self.setup,
            
            './zd/espresso/01-Density/scf.in': self.espresso_scf,
            
            './zd/espresso/02-Wfn/wfn.in': self.espresso_wfn,
            './zd/espresso/02-Wfn/wfn.pp.in': self.espresso_wfn_pp,
            './zd/espresso/02-Wfn/parabands.inp': self.espresso_wfn_parabands,

            './zd/espresso/03-Wfnq/wfn.in': self.espresso_wfnq,
            './zd/espresso/03-Wfnq/wfn.pp.in': self.espresso_wfnq_pp,

            './zd/espresso/04-Wfn_co/wfn.in': self.espresso_wfnco,
            './zd/espresso/04-Wfn_co/wfn.pp.in': self.espresso_wfnco_pp,
            './zd/espresso/04-Wfn_co/parabands.inp': self.espresso_wfnco_parabands,

            './zd/espresso/05-Wfn_fi/wfn.in': self.espresso_wfnfi,
            './zd/espresso/05-Wfn_fi/wfn.pp.in': self.espresso_wfnfi_pp,

            './zd/espresso/06-Wfnq_fi/wfn.in': self.espresso_wfnqfi,
            './zd/espresso/06-Wfnq_fi/wfn.pp.in': self.espresso_wfnqfi_pp,

            './zd/11-epsilon/epsilon.inp': self.epsilon,
            './zd/12-sigma/sigma.inp': self.sigma,
            './zd/13-kernel/kernel.inp': self.kernel,
            './zd/14-absorption/absorption.inp': self.absorption,

            './zd/finiteQ_grid/seed/13-kernel/kernel.inp': self.seed_kernel,
            './zd/finiteQ_grid/seed/14-absorption/absorption.inp': self.seed_absorption,

            './zd/epw/scf.in': self.epw_scf,
            './zd/epw/ph.in': self.epw_ph,
            './zd/epw/epw1.in': self.epw_xctph,
            './zd/epw/epw2.in': self.epw_ste,
            './zd/epw/epw3.in': self.epw_ste_eh_densities,

            './job_zd.sh': self.job_zd,
        }
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_zd.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './zd'
        ]

    def plot(self, **kwargs):
        EpwZdPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion