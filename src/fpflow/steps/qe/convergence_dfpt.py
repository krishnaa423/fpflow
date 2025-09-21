#region modules
from fpflow.steps.convergencestep_base import ConvergenceBaseStep
from fpflow.plots.convergence_dfpt import QeConvergenceDfptPlot
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class QeConvergenceDfptStep(ConvergenceBaseStep):
    def __init__(self, **kwargs):
        super().__init__(subdir1='qe', subdir2='dfpt', **kwargs)

    def plot(self, **kwargs):
        QeConvergenceDfptPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion