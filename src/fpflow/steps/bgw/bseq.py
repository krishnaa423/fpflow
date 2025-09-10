#region modules
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.io.read_write import str_2_f
from fpflow.structure.kpts import Kpts
import jmespath
from typing import List 
from fpflow.io.read_write import str_2_f
import os 
from fpflow.steps.step import Step 
from fpflow.inputs.grammars.namelist import NamelistGrammar
from fpflow.inputs.grammars.bgw import BgwGrammar
import jmespath
from fpflow.io.update import update_dict
from fpflow.io.logging import get_logger
from fpflow.schedulers.scheduler import Scheduler
from importlib.util import find_spec
from fpflow.structure.qe.qe_struct import QeStruct
from fpflow.structure.kpts import Kpts
import numpy as np
#endregion

#region variables
#endregion

#region functions
#endregion

#region classes
class BgwBseqStep(Step):
    @property
    def job_bseq(self):
        pass

    def get_kernel_strings(self, Qpt):
        qshift: list[int] = Qpt

        kerneldict: dict = {
            'exciton_Q_shift': f"2 {qshift[0]} {qshift[1]} {qshift[2]}",
            'use_symmetries_coarse_grid': '',
            'number_val_bands': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_cond_bands': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'use_wfn_hdf5': '',
            'dont_check_norms': '',
        }

        # Update if needed. 
        update_dict(kerneldict, jmespath.search('bse.kernel.args', self.inputdict))

        input_kernel = BgwGrammar().write(kerneldict)

        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bse.kernel.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf ../../epsmat.h5 ./
ln -sf ../../eps0mat.h5 ./
ln -sf ../../{jmespath.search('bse.absorption.wfnco_link', self.inputdict)} ./WFN_co.h5 
ln -sf ../../{jmespath.search('bse.absorption.wfnqco_link', self.inputdict)} ./WFNq_co.h5 
{scheduler.get_exec_prefix()}kernel.cplx.x &> kernel.inp.out
'''
        job_kernel = file_string
        
        return input_kernel, job_kernel

    def get_absorption_strings(self, Qpt):
        qshift: list[int] = Qpt
        pol_dir: list[int] = jmespath.search('bse.absorption.pol_dir[*]', self.inputdict)

        absorptiondict: dict = {
            'exciton_Q_shift': f"2 {qshift[0]} {qshift[1]} {qshift[2]}",
            'use_symmetries_coarse_grid': '',
            'use_symmetries_fine_grid': '',
            'use_symmetries_shifted_grid': '',
            'number_val_bands_coarse': jmespath.search('bse.absorption.val_bands', self.inputdict),
            'number_val_bands_fine': jmespath.search('bse.absorption.val_bands', self.inputdict) - 1,
            'number_cond_bands_coarse': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'number_cond_bands_fine': jmespath.search('bse.absorption.cond_bands', self.inputdict),
            'degeneracy_check_override': '',
            'diagonalization': '',
            # 'use_elpa': '',  
            'use_momentum': '',  
            # 'use_velocity': '',  
            'polarization': ' '.join(list(map(str, pol_dir))),
            'eqp_co_corrections': '',
            'dump_bse_hamiltonian': '',
            'use_wfn_hdf5': '',
            'energy_resolution': 0.1,
            'write_eigenvectors': jmespath.search('bse.absorption.nxct', self.inputdict),
            'dont_check_norms': '',
        }

        # Update if needed. 
        update_dict(absorptiondict, jmespath.search('bse.absorption.args', self.inputdict))

        input_absorption = BgwGrammar().write(absorptiondict)

        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bse.absorption.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf ../../epsmat.h5 ./
ln -sf ../../eps0mat.h5 ./
ln -sf ../../eqp1.dat eqp_co.dat 
#ln -sf ../../bsemat.h5 ./
ln -sf ../../{jmespath.search('bse.absorption.wfnco_link', self.inputdict)} ./WFN_co.h5 
ln -sf ../../{jmespath.search('bse.absorption.wfnqco_link', self.inputdict)} ./WFNq_co.h5 
ln -sf ../../{jmespath.search('bse.absorption.wfnfi_link', self.inputdict)} ./WFN_fi.h5 
ln -sf ../../{jmespath.search('bse.absorption.wfnqfi_link', self.inputdict)} ./WFNq_fi.h5 
{scheduler.get_exec_prefix()}absorption.cplx.x &> absorption.inp.out
mv bandstructure.dat bandstructure_absorption.dat
'''
        job_absorption = file_string
        
        return input_absorption, job_absorption

    def get_plotxct_strings(self, Qpt):
        plotxctdict: dict = {
            'hole_position': ' '.join(list(map(str, jmespath.search('bse.plotxct.hole_position', self.inputdict)))),
            'supercell_size': ' '.join(list(map(str, jmespath.search('bse.plotxct.supercell_size', self.inputdict)))),
            'use_symmetries_fine_grid': '',
            'use_symmetries_shifted_grid': '',
            'plot_spin': 1,
            'plot_state': jmespath.search('bse.plotxct.xct_state', self.inputdict),
            'use_wfn_hdf5': '',
        }

        # Update if needed. 
        update_dict(plotxctdict, jmespath.search('bse.plotxct.args', self.inputdict))

        input_plotxct = BgwGrammar().write(plotxctdict)

        scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bse.kernel.job_info')

        file_string = f'''#!/bin/bash
{scheduler.get_script_header()}

ln -sf ../../{jmespath.search('bse.absorption.wfnfi_link', self.inputdict)} ./WFN_fi.h5 
ln -sf ../../{jmespath.search('bse.absorption.wfnqfi_link', self.inputdict)} ./WFNq_fi.h5 
{scheduler.get_exec_prefix()}plotxct.cplx.x &> plotxct.inp.out 
volume.py ./scf.in espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
rm -rf *.a3Dr
'''
        job_plotxct = file_string
        
        return input_plotxct, job_plotxct

    def create_inputs_bseq(self):
        #TODO: Copy pasted. Can refactor.
        self.bseq_for_xctph_link_str: str = '\n\n\n\n'

        os.system('mkdir -p ./bseq')
        os.system('mkdir -p ./bseq_for_xctph')
        os.chdir('./bseq')

        #Qpts.
        Qpts: list = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('bseq.qgrid[0]', self.inputdict),
                jmespath.search('bseq.qgrid[1]', self.inputdict),
                jmespath.search('bseq.qgrid[2]', self.inputdict),
            ],
            is_reduced=False,
        ).bseq_qpts

        for Qpt_idx, Qpt in enumerate(Qpts):
            Qpt0 = f'{Qpt[0]:15.10f}'.strip()
            Qpt1 = f'{Qpt[1]:15.10f}'.strip()
            Qpt2 = f'{Qpt[2]:15.10f}'.strip()
            dir_name = f'Q_{Qpt0}_{Qpt1}_{Qpt2}'
            os.system(f'mkdir -p {dir_name}')
            self.bseq_for_xctph_link_str += f'ln -sf ../bseq/{dir_name} ./bseq_for_xctph/Q_{str(Qpt_idx).strip()}\n'
            os.chdir(f'./{dir_name}')
            
            inp_ker, job_ker = self.get_kernel_strings(Qpt)
            str_2_f(inp_ker, 'kernel.inp')
            str_2_f(job_ker, 'job_kernel.sh')

            inp_abs, job_abs = self.get_absorption_strings(Qpt)
            str_2_f(inp_abs, 'absorption.inp')
            str_2_f(job_abs, 'job_absorption.sh')


            inp_plotxct, job_plotxct = self.get_plotxct_strings(Qpt)
            str_2_f(inp_plotxct, 'plotxct.inp')
            str_2_f(job_plotxct, 'job_plotxct.sh')
            os.system('chmod u+x ./*.sh')

            os.chdir('../')

        os.chdir('../')

    def create_job_bseq(self):
        '''
        Idea is to create a list with start and stop indices to control execution.
        '''
        bseq_scheduler: Scheduler = Scheduler.from_jmespath(self.inputdict, 'bseq.job_info')
        kernel_scheduler: Scheduler  = Scheduler.from_jmespath(self.inputdict, 'bse.kernel.job_info')
        absorption_scheduler: Scheduler  = Scheduler.from_jmespath(self.inputdict, 'bse.absorption.job_info')
        plotxct_scheduler: Scheduler  = Scheduler.from_jmespath(self.inputdict, 'bse.plotxct.job_info')

        Qpts = Kpts.from_kgrid(
            kgrid = [
                jmespath.search('bseq.qgrid[0]', self.inputdict),
                jmespath.search('bseq.qgrid[1]', self.inputdict),
                jmespath.search('bseq.qgrid[2]', self.inputdict),
            ],
            is_reduced=False,
        ).bseq_qpts
        Qpts = np.array(Qpts)
        job_bseq = '#!/bin/bash\n'
        job_bseq += f'{bseq_scheduler.get_script_header()}\n'

        job_bseq += "start=0\n"
        job_bseq += f"stop={Qpts.shape[0]}\n\n"
        job_bseq += f"size={Qpts.shape[0]}\n\n"

        # Create the list.
        job_bseq += 'folders=('
        for Qpt in Qpts:
            Qpt0 = f'{Qpt[0]:15.10f}'.strip()
            Qpt1 = f'{Qpt[1]:15.10f}'.strip()
            Qpt2 = f'{Qpt[2]:15.10f}'.strip()
            subdir_name = f'Q_{Qpt0}_{Qpt1}_{Qpt2}'
            dir_name = f'"./bseq/{subdir_name}" '
            job_bseq += dir_name
        job_bseq += ')\n\n'


        kernel_commands = \
f'''    ln -sf ../../epsmat.h5 ./
    ln -sf ../../eps0mat.h5 ./
    ln -sf ../../{jmespath.search('bse.absorption.wfnco_link', self.inputdict)} ./WFN_co.h5 
    ln -sf ../../{jmespath.search('bse.absorption.wfnqco_link', self.inputdict)} ./WFNq_co.h5 
    {kernel_scheduler.get_exec_prefix()}kernel.cplx.x &> kernel.inp.out
'''
        
        absorption_commands = \
f'''    ln -sf ../../epsmat.h5 ./
    ln -sf ../../eps0mat.h5 ./
    ln -sf ../../eqp1.dat eqp_co.dat 
    ln -sf ../../{jmespath.search('bse.absorption.wfnco_link', self.inputdict)} ./WFN_co.h5 
    ln -sf ../../{jmespath.search('bse.absorption.wfnqco_link', self.inputdict)} ./WFNq_co.h5 
    ln -sf ../../{jmespath.search('bse.absorption.wfnfi_link', self.inputdict)} ./WFN_fi.h5 
    ln -sf ../../{jmespath.search('bse.absorption.wfnqfi_link', self.inputdict)} ./WFNq_fi.h5 
    {absorption_scheduler.get_exec_prefix()}absorption.cplx.x &> absorption.inp.out
    mv bandstructure.dat bandstructure_absorption.dat
'''
        
        plotxct_commands = \
f'''    ln -sf ../../{jmespath.search('bse.absorption.wfnfi_link', self.inputdict)} ./WFN_fi.h5 
ln -sf ../../{jmespath.search('bse.absorption.wfnqfi_link', self.inputdict)} ./WFNq_fi.h5 
    {plotxct_scheduler.get_exec_prefix()}plotxct.cplx.x &> plotxct.inp.out 
    volume.py ../../scf.in espresso *.a3Dr a3dr plotxct_elec.xsf xsf false abs2 true 
    rm -rf *.a3Dr
'''
        
        folder_variable = '${folders[$i]}'
        kpt_variable = '${i}'

        # Add the looping block.
        job_bseq += \
f'''
rm -rf ./bseq.out
touch ./bseq.out

LOG_FILE="$(pwd)/bseq.out"
exec &> "$LOG_FILE"

{self.bseq_for_xctph_link_str}

for (( i=$start; i<$stop; i++ )); do
    cd {folder_variable}

    echo -e "\\n\\n\\n"

    echo "Running {kpt_variable} th kpoint"
    echo "Entering folder {folder_variable}"
    
    echo "Starting kernel for {folder_variable}"
{kernel_commands}
    echo "Done kernel for {folder_variable}"

    echo "Starting absorption for {folder_variable}"
{absorption_commands}
    echo "Done absorption for {folder_variable}"

    echo "Starting plotxct for {folder_variable}"
{plotxct_commands}
    echo "Done plotxct for {folder_variable}"
    cd ../../

    echo "Exiting folder {folder_variable}"
done
'''

        str_2_f(job_bseq, 'job_bseq.sh')

    def create(self):
        self.create_inputs_bseq()
        self.create_job_bseq()
        os.system('chmod u+x ./*.sh')
    
    @property
    def job_scripts(self) -> List[str]:
        return [
            './job_bseq.sh'
        ]

    @property
    def save_inodes(self) -> List[str]:
        return []
    
    @property
    def remove_inodes(self) -> List[str]:
        return [
            './bseq',
            './bseq.out',
            './bseq_for_xctph',
            './job_bseq.sh',
        ]

#endregion