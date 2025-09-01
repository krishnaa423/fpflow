#region modules
import jmespath
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class JobInfo:
    def __init__(
        self,
        nodes: int = None,
        ntasks: int = None,
        time: int = None,
        nk: int = None,
        ni: int = None,
        **kwargs,
    ):
        self.nodes: int = nodes
        self.ntasks: int = ntasks
        self.time: str = time
        self.nk: int = nk
        self.ni: int = ni

        for key, value in kwargs.items():
            setattr(self, key, value)

    @staticmethod
    def from_inputdict(jmespath_str: str, inputdict: dict):
        info = jmespath.search(jmespath_str, inputdict)

        assert info is not None, f'The jmespath {jmespath_str} does not exist'

        return JobInfo(**info)
    
#endregion