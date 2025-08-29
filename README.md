# First-principles workflow

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
- ml_deepmd
- ml_dft
- ml_gw
- ml_bse
- create_script
- run_script
- remove_script
- plot_script
- interactive_script

## Steps for adding a new step.
- Fill out the step map in fpflow.steps.steps_map file. 
- Write the subclass of Step in fpflow.steps folder.

## How to work with the cli script. 
- get template input yaml: `fpflow --input=input_Si.yaml`. Or any other stored templates. 
- generate: `fpflow --generator=create`
- remove: `fpflow --generator=remove`
- run steps: `fpflow --manager=interactive|background`