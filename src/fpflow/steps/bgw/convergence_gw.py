#region modules
from fpflow.steps.convergencestep_base import ConvergenceBaseStep
from fpflow.plots.convergence_gw import BgwConvergenceGwPlot
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class BgwConvergenceGwStep(ConvergenceBaseStep):
    def __init__(self, **kwargs):
        super().__init__(subdir1='qe', subdir2='gw', **kwargs)

    def plot(self, **kwargs):
        BgwConvergenceGwPlot(inputdict=self.inputdict).save_figures(**kwargs)

#endregion