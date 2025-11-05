#region modules
from ase import Atoms 
import subprocess
import numpy as np 
from fpflow.io.read_write import str_2_f
import os 
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
def get_ibz_kpts(atoms: Atoms, kgrid: list, qshift: list = [0.0, 0.0, 0.0]):
    cell: np.ndarray = atoms.cell.array
    positions: np.ndarray = atoms.get_positions()
    numbers_idxs: list = []
    _, numbers_idxs = np.unique(atoms.numbers, return_inverse=True)
    numbers_idxs = (numbers_idxs + 1).tolist()

    positions_str: str = ''
    counter = 0
    for numbers_idx, position in zip(numbers_idxs, positions):
        counter += 1
        positions_str += f'{numbers_idx} {position[0]:15.10f} {position[1]:15.10f} {position[2]:15.10f}'
        if counter!=len(positions):
            positions_str += '\n'
    

    kgridinp_str: str = f'''{kgrid[0]} {kgrid[1]} {kgrid[2]}
0.0 0.0 0.0
{qshift[0]} {qshift[1]} {qshift[2]}
{cell[0][0]:15.10f} {cell[0][1]:15.10f} {cell[0][2]:15.10f}
{cell[1][0]:15.10f} {cell[1][1]:15.10f} {cell[1][2]:15.10f}
{cell[2][0]:15.10f} {cell[2][1]:15.10f} {cell[2][2]:15.10f}
{len(atoms)}
{positions_str}
0 0 0
.false.
.false.
.false.
'''
    
    str_2_f(kgridinp_str, 'kgrid.inp')
    result = subprocess.run([
        'kgrid.x',
        'kgrid.inp',
        'kgrid.log',
        'kgrid.out',
    ])

    kpts: np.ndarray = np.loadtxt('kgrid.log', skiprows=2)

    # Save the input and output files in a kgrid directory. 
    kgrid_folder_name = f'kgrid_{kgrid[0]}_{kgrid[1]}_{kgrid[2]}_qshift_{qshift[0]}_{qshift[1]}_{qshift[2]}'
    os.system(f'mkdir -p ./{kgrid_folder_name}')
    os.system(f'mv ./kgrid.inp ./{kgrid_folder_name}/kgrid.inp')
    os.system(f'mv ./kgrid.log ./{kgrid_folder_name}/kgrid.log')
    os.system(f'mv ./kgrid.out ./{kgrid_folder_name}/kgrid.out')

    return kpts
#endregion