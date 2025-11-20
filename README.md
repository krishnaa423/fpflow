# First-principles workflow

## Requirements

The `pip install -e .` will take care of installing the pip dependencies. Additionally, install 
Julia with `WannierIO.jl` package. This can be installed using 

```
using Pkg
Pkg.add("WannierIO")
```

This will be used for reading and writing data from wannier90 output. 

## Capabilities

Helps create input files, job scripts, manage runs, plots, and analysis. 

Supported software:
- Quantum Espresso
- BerkeleyGW
- Abacus
- Siesta
- Pyscf
- Gpaw

Steps, post-processsing, and (Plots):
- pseudos_qe
- esr_gen
- relax_qe
- cdft_qe
- md_qe
- scf_qe
- scf_abacus
- scf_siesta
- dftelbands_qe
- nscf_qe
- nscf_abacus
- nscf_siesta
- dos_qe 
- pdos_qe 
- kpdos_qe 
- wannier_qe
- dfpt_qe
- dfpt_elec_qe
- pp_dfpt_qe
- phbands_qe 
- phdos_qe 
- phmodes_qe  
- elph_grid_epw
- elph_bands_epw
- elph_frohlich_epw
- pp_elph_epw
- pol_epw
- pol_fp
- dmc_elph_fp
- phonopy_qe
- epsilon_bgw 
- sigmasc_bgw. Repeat this and below as many times as needed for scGW. 
- sigma_bgw
- gwelbands_bgw 
- kernel_bgw
- absorption_bgw 
- plotxct_bgw 
- bseq_bgw
- wannier_bgw
- xctph_epw 
- xctph_fp 
- xctpol_fp
- ste_epw
- zd_epw
- ste_fp  
- dmc_xctph_fp
- esf_fp
- esr_fp 
- ml_deepmd
- ml_dftqe
- ml_dfptqe
- ml_dvscqe
- ml_dnvscqe
- ml_dftabacus
- ml_dftsiesta
- ml_gwqe
- ml_bseqe
- convergence_scf_qe 
- convergence_dfpt_qe 
- convergence_gw_qe
- convergence_bse_qe 
- convergence_xctph
- convergence_xctpol
- convergence_ste
- convergence_dmc_elph
- convergence_dmc_xctph
- create_script
- run_script
- tail_script
- remove_script
- plot_script
- interactive_script

## Installation. 
Clone this repository. Then,

```
cd fpflow
pip install -e .
```


## Steps for adding a new step.
- Fill out the step map in fpflow.steps.steps_map file. 
- Write the subclass of Step in fpflow.steps folder.

## How to work with the cli script. 
- `fpflow input --list`: List all templates
- `fpflow input --template <template name>`: Generate input.yaml from template. 
- `fpflow generator --create`: Generate all the input files and job scripts.
- `fpflow generator --remove`: Remove input files and job scripts in the directory.
- `fpflow manager --run=interactive|background`: Run the job scripts in interactive mode or in the background. 
Only worth running it in the interactive queue either way.
- `fpflow manager --plot=no-gui|gui`: Plot data after runs and put them in the `./plots` subfolder. 