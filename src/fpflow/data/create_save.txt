#region: Modules.
from builtins import input
import numpy as np
import glob 
import os
#endregion

#region: Variables.
#endregion

#region: Functions.
def get_nqpt():
    fname = './tmp/_ph0/struct.phsave/control_ph.xml'

    fid = open(fname, 'r')
    lines = fid.readlines()
    # these files are relatively small so reading the whole thing shouldn't
    # be an issue
    fid.close()

    line_number_of_nqpt = 0
    while 'NUMBER_OF_Q_POINTS' not in lines[line_number_of_nqpt]:
        # increment to line of interest
        line_number_of_nqpt += 1
    line_number_of_nqpt += 1  # its on the next line after that text

    nqpt = int(lines[line_number_of_nqpt])

    return nqpt

def get_first_nonzero_dvscf_file(qpt_idx):

   if qpt_idx==1:
      return './tmp/_ph0/struct.dvscf1'
   else:
      glob_pat = f'./tmp/**/struct.q_{qpt_idx}/struct.dvscf1'

      filelist = glob.glob(glob_pat)

      for filename in filelist:
         filesize = os.path.getsize(filename)

         # Return if the filename is greater than zero. 
         if filesize > 0 : return filename 

def move_data_to_save():

   os.system('mkdir -p save')

   # Move the phsave folder. 
   os.system('cp -r ./tmp/_ph0/struct.phsave save/')

   # Move files related to each qpt. 
   for qpt_idx in range(1, get_nqpt()+1):

      # Move the dvscf file. 
      dvscf_file = get_first_nonzero_dvscf_file(qpt_idx)
      os.system(f'cp {dvscf_file} ./save/struct.dvscf_q{qpt_idx}')

      # Move the dyn files. 
      os.system(f'cp struct.dyn{qpt_idx} save/struct.dyn_q{qpt_idx}')

      # Delete some files. 
      # os.system(f'rm ./tmp/**/struct.q_{qpt_idx}/*wfc*' )

def main():
    move_data_to_save()
#endregion

#region: Classes.
#endregion

#region: Main.
if __name__=='__main__':
    main()
#endregion

#region: Old
# prefix='struct'
# nqpt=get_nqpt(prefix)

# os.system('mkdir save')

# for iqpt in np.arange(1,nqpt+1):
#   label = str(iqpt)

#   os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
#   if (iqpt == 1):
#     os.system('cp ./tmp/_ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
#     os.system('cp -r ./tmp/_ph0/'+prefix+'.phsave save/')
#   else:
#     os.system('cp ./tmp/_ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
#     os.system('rm ./tmp/_ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
#endregion: old