#region modules
from fpflow.steps.convergencestep_base import ConvergenceBaseStep
from fpflow.plots.convergence_scf import QeConvergenceScfPlot

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeConvergenceScfStep(ConvergenceBaseStep):
    def __init__(self, **kwargs):
        super().__init__(subdir1='qe', subdir2='scf', **kwargs)

    def plot(self, **kwargs):
        QeConvergenceScfPlot().save_figures(**kwargs)

#endregion