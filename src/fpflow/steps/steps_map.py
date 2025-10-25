#region modules
from fpflow.steps.qe.pseudos import QePseudosStep
from fpflow.steps.esr_gen import EsrgenStep
from fpflow.steps.qe.relax import QeRelaxStep
from fpflow.steps.qe.cdft import QeCdftStep
from fpflow.steps.qe.md import QeMdStep
from fpflow.steps.qe.scf import QeScfStep
from fpflow.steps.qe.scf_pw2bgw import QeScfPw2bgwStep
from fpflow.steps.qe.dfpt import QeDfptStep
from fpflow.steps.qe.phbands import QePhbandsStep
from fpflow.steps.qe.phdos import QePhdosStep
from fpflow.steps.qe.phmodes import QePhmodesStep
from fpflow.steps.qe.dos import QeDosStep
from fpflow.steps.qe.pdos import QePdosStep
from fpflow.steps.qe.dftelbands import QeDftelbandsStep
from fpflow.steps.qe.kpdos import QeKpdosStep
from fpflow.steps.qe.wfnwannier import QeWfnWannierStep
from fpflow.steps.qe.wannier import QeWannierStep
from fpflow.steps.qe.wfn import QeWfnStep
from fpflow.steps.qe.epw import QeEpwStep
from fpflow.steps.qe.elph import QeElphEpwStep
from fpflow.steps.qe.wfnq import QeWfnqStep
from fpflow.steps.qe.wfnfi import QeWfnfiStep
from fpflow.steps.qe.wfnqfi import QeWfnqfiStep
from fpflow.steps.qe.phonopy import QePhonopyStep
from fpflow.steps.qe.phonopy_nested import QePhonopyNestedStep
from fpflow.steps.bgw.epsilon import BgwEpsilonStep
from fpflow.steps.bgw.sigma import BgwSigmaStep
from fpflow.steps.bgw.gwelbands import BgwGwelbandsStep
from fpflow.steps.bgw.kernel import BgwKernelStep
from fpflow.steps.bgw.absorption import BgwAbsorptionStep
from fpflow.steps.bgw.plotxct import BgwPlotxctStep
from fpflow.steps.bgw.bseq import BgwBseqStep
from fpflow.steps.bgw.bseq_epw import EpwBseqStep
from fpflow.steps.xctph import XctphStep
from fpflow.steps.bgw.xctph_epw import EpwXctphStep
from fpflow.steps.xctpol import XctpolStep
from fpflow.steps.bgw.ste_epw import EpwSteStep
from fpflow.steps.ste import SteStep
from fpflow.steps.qe.zd_epw import EpwZdStep
from fpflow.steps.esf import EsfStep
from fpflow.steps.esr import EsrStep
from fpflow.steps.qe.convergence_scf import QeConvergenceScfStep
from fpflow.steps.qe.convergence_dfpt import QeConvergenceDfptStep
from fpflow.steps.bgw.convergence_gw import BgwConvergenceGwStep
from fpflow.steps.bgw.convergence_bse import BgwConvergenceBseStep
from fpflow.steps.ml.deepmd import MlDeepmdStep
from fpflow.steps.ml.qe.dft import MlQeDftStep
from fpflow.steps.ml.qe.gw import MlQeGwStep
from fpflow.steps.ml.qe.bse import MlQeBseStep
from fpflow.steps.post_steps.create_script import CreateScriptStep
from fpflow.steps.post_steps.run_script import RunScriptStep
from fpflow.steps.post_steps.tail_script import TailScriptStep
from fpflow.steps.post_steps.remove_script import RemoveScriptStep
from fpflow.steps.post_steps.plot_script import PlotScriptStep
from fpflow.steps.post_steps.interactive_script import InteractiveScriptStep
#endregion

#region variables
step_class_map: dict = {
    'pseudos_qe': QePseudosStep,
    'esr_gen': EsrgenStep,
    'relax_qe': QeRelaxStep,
    'cdft_qe': QeCdftStep,
    'md_qe': QeMdStep,
    'ml_deepmd': MlDeepmdStep, #TODO
    'scf_qe': QeScfStep,
    'scf_pw2bgw_qe': QeScfPw2bgwStep,
    'scf_abacus': QeScfStep, #TODO
    'scf_siesta': QeScfStep, #TODO
    'ml_dftqe': MlQeDftStep,
    'ml_dfptqe': MlQeDftStep, #TODO
    'ml_dvscqe': MlQeDftStep, #TODO
    'ml_dnvscqe': MlQeDftStep, #TODO
    'ml_dftabacus': MlQeDftStep, #TODO
    'ml_dftsiesta': MlQeDftStep, #TODO
    'dfpt_qe': QeDfptStep,
    'pp_dvscf_dfpt': QeDfptStep, #TODO
    'elph_dfpt': QeDfptStep, #TODO
    'phbands_qe': QePhbandsStep,
    'phdos_qe': QePhdosStep, #TODO
    'phmodes_qe': QePhmodesStep,
    'dos_qe': QeDosStep, 
    'pdos_qe': QePdosStep, 
    'dftelbands_qe': QeDftelbandsStep,
    'kpdos_qe': QeKpdosStep, 
    'wfnwannier_qe': QeWfnWannierStep, 
    'wannier_qe': QeWannierStep,
    'wfn_qe': QeWfnStep,
    'pp_wfnqe': QeWfnStep, #TODO
    'pp_wfnqesym': QeWfnStep, #TODO
    'wfn_abacus': QeWfnStep, #TODO
    'wfn_siesta': QeWfnStep, #TODO
    'epw_qe': QeEpwStep, 
    'elph_epw': QeElphEpwStep, #TODO
    'pol_fp': QeElphEpwStep, #TODO
    'wfnq_qe': QeWfnqStep,
    'wfnq_abacus': QeWfnqStep, #TODO
    'wfnq_siesta': QeWfnqStep, #TODO
    'wfnfi_qe': QeWfnfiStep,
    'wfnqfi_qe': QeWfnqfiStep,
    'phonopy_qe': QePhonopyNestedStep,
    'pp_dvscf_phonopy': QePhonopyStep, #TODO
    'elph_phonopy': QePhonopyStep, #TODO
    'epsilon_bgw': BgwEpsilonStep,
    'sigmasc_bgw': BgwSigmaStep, #TODO
    'sigma_bgw': BgwSigmaStep, 
    'ml_gwqe': MlQeGwStep, #TODO
    'gwelbands_bgw': BgwGwelbandsStep,
    'kernel_bgw': BgwKernelStep, 
    'absorption_bgw': BgwAbsorptionStep,
    'plotxct_bgw': BgwPlotxctStep,
    'bseq_bgw': BgwBseqStep,
    'bseq_epw': EpwBseqStep,
    'ml_bseqe': MlQeBseStep, #TODO
    'xctwannier_bgw': MlQeBseStep, #TODO
    'xctph_fp': XctphStep, #TODO
    'xctph_epw': EpwXctphStep,
    'xctpol_fp': XctpolStep, #TODO
    'ste_epw': EpwSteStep, 
    'ste_fp': SteStep, #TODO
    'zd_epw': EpwZdStep,
    'esf_fp': EsfStep, #TODO
    'esr_fp': EsrStep,  #TODO
    'convergence_scf_qe': QeConvergenceScfStep, 
    'convergence_dfpt_qe': QeConvergenceDfptStep, 
    'convergence_gw_qe': BgwConvergenceGwStep, 
    'convergence_bse_qe': BgwConvergenceBseStep, 
    'create_script': CreateScriptStep,
    'run_script': RunScriptStep,
    'tail_script': TailScriptStep,
    'remove_script': RemoveScriptStep,
    'plot_script': PlotScriptStep,
    'interactive_script': InteractiveScriptStep,
}
#endregion

#region functions
#endregion

#region classes
#endregion