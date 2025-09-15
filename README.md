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
- scf_abacus
- scf_siesta
- ml_dftqe
- ml_dfptqe
- ml_dvscqe
- ml_dnvscqe
- ml_dftabacus
- ml_dftsiesta
- dfpt_qe
- pp_dvscf_dfpt
- elph_dfpt (plot: elph_proj) 
- phbands_qe (plot: phbands)
- phdos_qe (plot: phdos)
- phmodes_qe (plot: phmode)
- dos_qe (plot: dos)
- pdos_qe (plot: dos)
- dftelbands_qe (plot: dftelbands)
- kpdos_qe (plot: kpdos)
- wannier_qe (plot: wannier_bands, wannier_wfn, wannier_hr, wannier_hr_decay)
- tbmodels
- wfn_qe.  Use qe, abacus or siesta. (plot: wfn after post processing) 
- pp_wfnqe
- pp_wfnqesym
- wfn_abacus
- wfn_siesta
- epw_qe
- elph_epw (plot: elph_proj) 
- wfnq_qe.  Use qe, abacus or siesta. 
- wfnq_abacus
- wfnq_siesta
- wfnfi_qe
- wfnqfi_qe
- phonopy_qe (plot: phonopy_bands)
- pp_dvscf_phonopy
- elph_phonopy (plot: elph_proj) 
- phonopy_abacus
- phonopy_siesta
- epsilon_bgw (plot: conv, epsmat)
- sigmasc_bgw. Repeat this and below as many times as needed for scGW. 
- epsilonsc_bgw
- sigma_bgw (plot: conv, sigma)
- ml_gwqe
- gwelbands_bgw (plot: gwbands, dftelbands)
- kernel_bgw
- absorption_bgw (plot: absorption, hbse)
- plotxct_bgw (plot: xct_wfn)
- bseq_bgw (plot: xctbands)
- ml_bseqe
- xctwannier_bgw
- xctph (plot: xctph_proj)
- xctpol (plot: xct_shifts, ph_shifts, disp, ph_proj, wfn, conv, emission)
- ste (plot: xct_shifts, ph_shifts, disp, ph_proj, wfn, conv, emission, barrier)
- esf (plot: xsf with forces, ph_proj)
- esr (plot: e_change, *_change)
- convergence_scf_qe (plot: dftelbands_ecut_kgrid)
- convergence_dfpt_qe (plot: phbands_qgrid_tr)
- convergence_gw_qe (plot: gwelbands_ecut_kgrid_bands)
- convergence_bse_qe (plot: absorption_cvk)
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