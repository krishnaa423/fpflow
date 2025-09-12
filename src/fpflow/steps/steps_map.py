#region modules
from fpflow.steps.qe.pseudos import QePseudosStep
from fpflow.steps.esr_gen import EsrgenStep
from fpflow.steps.qe.relax import QeRelaxStep
from fpflow.steps.qe.cdft import QeCdftStep
from fpflow.steps.qe.md import QeMdStep
from fpflow.steps.qe.scf import QeScfStep
from fpflow.steps.qe.dfpt import QeDfptStep
from fpflow.steps.qe.phbands import QePhbandsStep
from fpflow.steps.qe.phdos import QePhdosStep
from fpflow.steps.qe.phmodes import QePhmodesStep
from fpflow.steps.qe.dos import QeDosStep
from fpflow.steps.qe.pdos import QePdosStep
from fpflow.steps.qe.dftelbands import QeDftelbandsStep
from fpflow.steps.qe.kpdos import QeKpdosStep
from fpflow.steps.qe.wannier import QeWannierStep
from fpflow.steps.qe.wfn import QeWfnStep
from fpflow.steps.qe.epw import QeEpwStep
from fpflow.steps.qe.elph import QeElphStep
from fpflow.steps.qe.wfnq import QeWfnqStep
from fpflow.steps.qe.wfnfi import QeWfnfiStep
from fpflow.steps.qe.wfnqfi import QeWfnqfiStep
from fpflow.steps.qe.phonopy import QePhonopyStep
from fpflow.steps.bgw.epsilon import BgwEpsilonStep
from fpflow.steps.bgw.sigma import BgwSigmaStep
from fpflow.steps.bgw.gwelbands import BgwGwelbandsStep
from fpflow.steps.bgw.kernel import BgwKernelStep
from fpflow.steps.bgw.absorption import BgwAbsorptionStep
from fpflow.steps.bgw.plotxct import BgwPlotxctStep
from fpflow.steps.bgw.bseq import BgwBseqStep
from fpflow.steps.xctph import XctphStep
from fpflow.steps.xctpol import XctpolStep
from fpflow.steps.ste import SteStep
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
    'scf_qe': QeScfStep,
    'dfpt_qe': QeDfptStep,
    'phbands_qe': QePhbandsStep,
    'phdos_qe': QePhdosStep,
    'phmodes_qe': QePhmodesStep,
    'dos_qe': QeDosStep,
    'pdos_qe': QePdosStep,
    'dftelbands_qe': QeDftelbandsStep,
    'kpdos_qe': QeKpdosStep,
    'wannier_qe': QeWannierStep,
    'wfn_qe': QeWfnStep,
    'epw_qe': QeEpwStep,
    'elph': QeElphStep,
    'wfnq_qe': QeWfnqStep,
    'wfnfi_qe': QeWfnfiStep,
    'wfnqfi_qe': QeWfnqfiStep,
    'phonopy_qe': QePhonopyStep,
    'epsilon_bgw': BgwEpsilonStep,
    'sigma_bgw': BgwSigmaStep,
    'gwelbands_bgw': BgwGwelbandsStep,
    'kernel_bgw': BgwKernelStep,
    'absorption_bgw': BgwAbsorptionStep,
    'plotxct_bgw': BgwPlotxctStep,
    'bseq_bgw': BgwBseqStep,
    'xctph': XctphStep,
    'xctpol': XctpolStep,
    'ste': SteStep,
    'esf': EsfStep,
    'esr': EsrStep,
    'convergence_scf_qe': QeConvergenceScfStep,
    'convergence_dfpt_qe': QeConvergenceDfptStep,
    'convergence_gw_qe': BgwConvergenceGwStep,
    'convergence_bse_qe': BgwConvergenceBseStep,
    'mldeepmd': MlDeepmdStep,
    'mldftqe': MlQeDftStep,
    'mlgwqe': MlQeGwStep,
    'mlbseqe': MlQeBseStep,
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