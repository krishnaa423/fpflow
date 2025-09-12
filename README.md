# First-principles workflow

Helps create input files, job scripts, manage runs, plots, and analysis. 

Supported software:
- Quantum Espresso
- BerkeleyGW
- Abacus
- Siesta
- Pyscf
- Gpaw

Steps:
- pseudos_qe
- esr_gen
- relax_qe
- cdft_qe
- md_qe
- scf_qe
- dfpt_qe
- phbands_qe
- phdos_qe
- phmodes_qe
- dos_qe
- pdos_qe
- dftelbands_qe
- kpdos_qe
- wannier_qe
- wfn_qe
- epw_qe
- elph
- wfnq_qe
- wfnfi_qe
- wfnqfi_qe
- phonopy_qe
- epsilon_bgw
- sigma_bgw
- kernel_bgw
- absorption_bgw
- plotxct_bgw
- bseq_bgw
- xctph
- xctpol
- ste
- esf
- convergence_scf
- convergence_dfpt
- convergence_gw
- convergence_bse
- ml_deepmd
- ml_dft
- ml_gw
- ml_bse
- create_script
- run_script
- remove_script
- plot_script
- interactive_script

## Installation. 

```
git clone 
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