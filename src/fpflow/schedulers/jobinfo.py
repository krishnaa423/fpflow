#region modules
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
    def from_job_id(id: str, input_dict: dict):
        assert id in input_dict['job_types'], f'id: {id} is not in job_types'

        return JobInfo(**input_dict['job_types'][id])
    
#endregion