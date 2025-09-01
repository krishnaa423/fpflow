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
from fpflow.steps.bgw.kernel import BgwKernelStep
from fpflow.steps.bgw.absorption import BgwAbsorptionStep
from fpflow.steps.bgw.plotxct import BgwPlotxctStep
from fpflow.steps.bgw.bseq import BgwBseqStep
from fpflow.steps.xctph import XctphStep
from fpflow.steps.xctpol import XctpolStep
from fpflow.steps.ste import SteStep
from fpflow.steps.esf import EsfStep
from fpflow.steps.ml.deepmd import MlDeepmdStep
from fpflow.steps.ml.dft import MlDftStep
from fpflow.steps.ml.gw import MlGwStep
from fpflow.steps.ml.bse import MlBseStep
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
    'kernel_bgw': BgwKernelStep,
    'absorption_bgw': BgwAbsorptionStep,
    'plotxct_bgw': BgwPlotxctStep,
    'bseq_bgw': BgwBseqStep,
    'xctph': XctphStep,
    'xctpol': XctpolStep,
    'ste': SteStep,
    'esf': EsfStep,
    'ml_deepmd': MlDeepmdStep,
    'ml_dft': MlDftStep,
    'ml_gw': MlGwStep,
    'ml_bse': MlBseStep,
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