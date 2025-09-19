#region modules
from fpflow.steps.convergencestep_base import ConvergenceBaseStep
from fpflow.plots.convergence_bse import BgwConvergenceBsePlot
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class BgwConvergenceBseStep(ConvergenceBaseStep):
    def __init__(self, **kwargs):
        super().__init__(subdir1='qe', subdir2='bse', **kwargs)

    def plot(self, **kwargs):
        BgwConvergenceBsePlot().save_figures(**kwargs)

#endregion