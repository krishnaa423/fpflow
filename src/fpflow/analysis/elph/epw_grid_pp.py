#region modules
import h5py 

#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class EpwElphGridPp:
    '''
    run() is executed in the ./elph_grid_epw directory. 
    '''
    def __init__(self):
        self.datasets: dict = {}

    def read_coarse_elph(self):
        pass

    def read_coarse_dynq(self):
        pass

    def read_coarse_pheigs(self):
        pass

    def read_coarse_phevecs(self):
        pass

    def read_coarse_born(self):
        pass

    def read_coarse_eps(self):
        pass

    def read_coarse(self):
        self.read_coarse_elph()
        self.read_coarse_dynq()
        self.read_coarse_pheigs()
        self.read_coarse_phevecs()
        self.read_coarse_born()
        self.read_coarse_eps()

    def read_fine_elph(self):
        pass

    def read_fine_dynq(self):
        pass

    def read_fine_pheigs(self):
        pass

    def read_fine_phevecs(self):
        pass

    def read_fine(self):
        self.read_fine_elph()
        self.read_fine_pheigs()
        self.read_fine_phevecs()
        self.read_fine_dynq()

    def read_all(self):
        self.read_coarse()
        self.read_fine()

    def write_all(self):
        with h5py.File('./elph.h5', 'w') as f:
            for name, data in self.datasets.items():
                f.create_dataset(name, data=data) 

    def run(self):
        self.read_all()
        self.write_all()

#endregion
