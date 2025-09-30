# First-principles workflow

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
- ml_deepmd
- scf_qe. Use qe, abacus or siesta. 
- scf_pw2bgw_qe
- scf_abacus
- scf_siesta
- ml_dftqe
- ml_dfptqe
- ml_dvscqe
- ml_dnvscqe
- ml_dftabacus
- ml_dftsiesta
- dfpt_qe
- phbands_qe 
- phdos_qe 
- phmodes_qe 
- dos_qe 
- pdos_qe 
- dftelbands_qe 
- kpdos_qe 
- wannier_qe 
- tbmodels
- wfn_qe
- wfn_abacus
- wfn_siesta
- epw_qe
- elph_epw 
- pol_fp
- wfnq_qe
- wfnq_abacus
- wfnq_siesta
- wfnfi_qe
- wfnqfi_qe
- phonopy_qe 
- epsilon_bgw 
- sigmasc_bgw. Repeat this and below as many times as needed for scGW. 
- sigma_bgw
- ml_gwqe
- gwelbands_bgw 
- kernel_bgw
- absorption_bgw 
- bse_waninterp_bgw
- plotxct_bgw 
- bse_fp 
- bseq_bgw 
- bseq_waninterp_bgw
- bseq_fp 
- ml_bseqe
- wannier_bgw
- xctph_fp 
- xctpol_fp 
- xctpol_dmc_fp
- ste_fp 
- esf_fp
- esr_fp 
- convergence_scf_qe 
- convergence_dfpt_qe 
- convergence_gw_qe
- convergence_bse_qe 
- create_script
- run_script
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