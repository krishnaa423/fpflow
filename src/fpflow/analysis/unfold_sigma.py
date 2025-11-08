# File copied and modified from BerkeleyGW. 
#region modules
import numpy as np
import h5py

#endregion

#region variables
#endregion

#region functions
def unfold_sigma():
  ''' 
  Unfold sigma_hp.log and eqp0.dat and eqp1.dat files using symmetries in kgrid.out
  Only nspin = 1 considered
  kpoints in sigma_hp.log must be identical (values and order) to the irreducible k-points in kgrid.log

  scissors operation (simplest linear version) provided for eqp.dat
  '''

  scissors = True  # linear scissors operation
  rydberg2eV = 13.60569253
  _tiny = 1E-8
  wfn_has_symmetry = True  # if WFN_inner has no symmetyr and therefore contains full k-points. 

  # process kgrid.log file

  fklog = open('./kgrid/kgrid.out', 'r')
  klines = fklog.readlines()
  fklog.close()

  for iline in range(len(klines)):
    ln = klines[iline].split()

    if len(ln) > 3 and ln[3] == 'original':
      nextln = klines[iline+1].split()
      nktot = int(nextln[0])

      ktot_ln_start = iline + 2  # record the line number to read full k-point list

    if len(ln) > 3 and ln[3] == 'irreducible':
      nextln = klines[iline+1].split()
      nksym = int(nextln[0])

      ksym_ln_start = iline + 2  # record the line number to read full k-point list

  # now prepare and read the klists and the mapping
  ktot = np.zeros((nktot, 3), dtype='f8')
  ksym = np.zeros((nksym, 3), dtype='f8')
  ktotmap = np.zeros(nktot, dtype='i4')
  kmap_sym_to_tot = np.zeros(nksym, dtype='i4')
  kmap_tot_to_sym = np.zeros(nktot, dtype='i4')

  print('Number of k-points in full grid:', nktot)
  print('Number of k-points in irreducible wedge:', nksym)

  for ik in range(nktot):
    ln = klines[ik + ktot_ln_start].split()

    ktot[ik, 0] = float(ln[1])
    ktot[ik, 1] = float(ln[2])
    ktot[ik, 2] = float(ln[3])

    kidx = int(ln[5])

    if kidx == 0:
      ktotmap[ik] = ik
    else:
      ktotmap[ik] = kidx - 1

  for ik in range(nksym):
    ln = klines[ik + ksym_ln_start].split()

    ksym[ik, 0] = float(ln[1])
    ksym[ik, 1] = float(ln[2])
    ksym[ik, 2] = float(ln[3])

    kmap_sym_to_tot[ik] = int(ln[5]) - 1

  # generate needed kmap
  kmap_tot_to_sym[:] = -1
  for ik in range(nktot):
    for iksym in range(nksym):
      if ktotmap[ik] == kmap_sym_to_tot[iksym]:
        kmap_tot_to_sym[ik] = iksym

  # process sigma_hp.log file
  fsig = open('./sigma_hp.log', 'r')
  siglines = fsig.readlines()
  fsig.close()

  for iline in range(len(siglines)):
    ln = siglines[iline].split()

    if len(ln) > 0 and ln[0] == 'band_index':
      bmin = int(ln[1])
      bmax = int(ln[2])

      nband = bmax - bmin + 1

  print('Number of bands:', nband)

  data = []

  for iline in range(len(siglines)):
    ln = siglines[iline].split()

    if len(ln) > 0 and ln[0] == 'k':
      bands = []

      for ib in range(nband):
        bands.append(siglines[iline + 3 + ib])

      data.append(bands)

  if len(data) != nksym:
    raise ValueError('Number of kpoints in reduced BZ mismatch in kgrid.log and sigma_hp.log.')

  # output sigma_hp
  fsig_unfold = open('./sigma_hp_unfold.log', 'w')

  for iline in range(len(siglines)):
    ln = siglines[iline].split()
    if len(ln) > 0 and ln[0] == 'k':
      break

    fsig_unfold.write(siglines[iline])

  for ik in range(nktot):
    fsig_unfold.write('       k =  ')
    fsig_unfold.write('%.6f'%ktot[ik,0] + '  ' + '%.6f'%ktot[ik,1] + '  ' + '%.6f'%ktot[ik,2])
    fsig_unfold.write(' ik =   0 spin = 1\n\n')  # we do not trace ik index in WFN file
 
    fsig_unfold.write('   n            Emf             Eo              X           SX-X      ')
    fsig_unfold.write('       CH            Sig            Vxc           Eqp0        ')
    fsig_unfold.write('   Eqp1            Znk\n')
    for ib in range(nband):
      fsig_unfold.write(data[kmap_tot_to_sym[ik]][ib])

    fsig_unfold.write('\n')

  fsig_unfold.close()

  # if scissors, read DFT band energies from WFN_inner.h5
  if scissors:
    wfn = h5py.File('./WFN_inner.h5', 'r')

    print('reading band energies from WFN_inner.h5')

    bands_wfn = wfn['mf_header/kpoints/el'][()] * rydberg2eV
    kwfn = wfn['mf_header/kpoints/rk'][()]

    ns_wfn = len(bands_wfn)
    if ns_wfn > 1:
      raise ValueError('Only nspin=1 is allowed.')

    nk_wfn = len(bands_wfn[0])
    
    nb_wfn = len(bands_wfn[0,0])

    # now we generate kmap_tot_to_wfn
    kmap_tot_to_wfn = np.zeros(nktot, dtype='i4')

    # if WFN_inner contains only symmetry-reduced kpoints
    if wfn_has_symmetry:
      # we find the map between ksym to kwfn first
      kmap_sym_to_wfn = np.zeros(nksym, dtype='i4')
  
      for ik in range(nksym):
        for jk in range(nk_wfn):
          kdiff = ksym[ik] - kwfn[jk]
          if np.sqrt(np.dot(kdiff, kdiff)) < _tiny:
            kmap_sym_to_wfn[ik] = jk
            break
      
      for ik in range(nktot):
        kmap_tot_to_wfn[ik] = kmap_sym_to_wfn[kmap_tot_to_sym[ik]]

    # if WFN_inner contains all full k-grid
    else:
      for ik in range(nktot):
        for jk in range(nk_wfn):
          kdiff = ktot[ik] - kwfn[jk]
          if np.sqrt(np.dot(kdiff, kdiff)) < _tiny:
            kmap_tot_to_wfn[ik] = jk
            break

  # output eqp.dat
  feqp0 = open('./eqp0_unfold.dat', 'w')
  feqp1 = open('./eqp1_unfold.dat', 'w')

  # save gwcor for all kpoints
  gwcor0_lower = np.zeros(nktot, dtype=float)
  gwcor1_lower = np.zeros(nktot, dtype=float)
  gwcor0_upper = np.zeros(nktot, dtype=float)
  gwcor1_upper = np.zeros(nktot, dtype=float)

  for ik in range(nktot):
    feqp0.write('  ' + '%.8f'%ktot[ik,0] + '  ' + '%.8f'%ktot[ik,1] + '  ' + '%.9f'%ktot[ik,2] + '  ')
    feqp0.write('      ' + str(nband) + '\n')

    feqp1.write('  ' + '%.8f'%ktot[ik,0] + '  ' + '%.8f'%ktot[ik,1] + '  ' + '%.9f'%ktot[ik,2] + '  ')
    feqp1.write('      ' + str(nband) + '\n')

    for ib in range(nband):
      ln = data[kmap_tot_to_sym[ik]][ib].split()
      edft = float(ln[1])
      eqp0 = float(ln[8])
      eqp1 = float(ln[9])

      # find gwcor for lower side and upper side
      # independent of whether using scissors
      if ib == 0:
        gwcor0_lower[ik] = eqp0 - edft
        gwcor1_lower[ik] = eqp1 - edft
        #print('GW corrections on lower end for eqp0 and eqp1', gwcor0_lower[ik], gwcor1_lower[ik])

      if ib == nband - 1:
        gwcor0_upper[ik] = eqp0 - edft
        gwcor1_upper[ik] = eqp1 - edft
        #print('GW corrections on upper end for eqp0 and eqp1', gwcor0_upper[ik], gwcor1_upper[ik])

      # if scissors, cross-check DFT band energies
      if scissors:
        edft_wfn = bands_wfn[0, kmap_tot_to_wfn[ik], ib + bmin - 1]

        if np.abs(edft_wfn - edft) > _tiny:
          print('sigma_hp energy', edft, 'WFN energy', edft_wfn)
          raise ValueError('Inconsistent energy between sigma_hp.log and WFN_inner.h5')

      feqp0.write('       1')
      feqp0.write('%8d'%(bmin + ib) + '%15.9f'%edft + '%15.9f'%eqp0 + '\n')
  
      feqp1.write('       1')
      feqp1.write('%8d'%(bmin + ib) + '%15.9f'%edft + '%15.9f'%eqp1 + '\n')

  feqp0.close()
  feqp1.close()

  # write eqp.dat for scissors case
  if scissors:
    feqp0 = open('./eqp0_unfold_scissors.dat', 'w')
    feqp1 = open('./eqp1_unfold_scissors.dat', 'w')

    fsig_scissors = open('./sigma_hp_unfold_scissors.log', 'w')
    for iline in range(len(siglines)):
      ln = siglines[iline].split()
      if len(ln) > 0 and ln[0] == 'k':
        break

      fsig_scissors.write(siglines[iline])

    for ik in range(nktot):
      feqp0.write('  ' + '%.8f'%ktot[ik,0] + '  ' + '%.8f'%ktot[ik,1] + '  ' + '%.9f'%ktot[ik,2] + '  ')
      feqp0.write('      ' + str(nb_wfn) + '\n')
  
      feqp1.write('  ' + '%.8f'%ktot[ik,0] + '  ' + '%.8f'%ktot[ik,1] + '  ' + '%.9f'%ktot[ik,2] + '  ')
      feqp1.write('      ' + str(nb_wfn) + '\n')

      fsig_scissors.write('       k =  ')
      fsig_scissors.write('%.6f'%ktot[ik,0] + '  ' + '%.6f'%ktot[ik,1] + '  ' + '%.6f'%ktot[ik,2])
      fsig_scissors.write(' ik =   0 spin = 1\n\n')  # we do not trace ik index in WFN file
   
      fsig_scissors.write('   n            Emf             Eo              X           SX-X      ')
      fsig_scissors.write('       CH            Sig            Vxc           Eqp0        ')
      fsig_scissors.write('   Eqp1            Znk\n')

      for ib_wfn in range(nb_wfn):
        edft = bands_wfn[0, kmap_tot_to_wfn[ik], ib_wfn]

        if ib_wfn < bmin-1:
          eqp0 = edft + gwcor0_lower[ik]
          eqp1 = edft + gwcor1_lower[ik]
        elif ib_wfn > bmax-1:
          eqp0 = edft + gwcor0_upper[ik]
          eqp1 = edft + gwcor1_upper[ik]
        else:
          ln = data[kmap_tot_to_sym[ik]][ib_wfn + 1 - bmin].split()
          edft = float(ln[1])
          eqp0 = float(ln[8])
          eqp1 = float(ln[9])

        feqp0.write('       1')
        feqp0.write('%8d'%(ib_wfn + 1) + '%15.9f'%edft + '%15.9f'%eqp0 + '\n')
  
        feqp1.write('       1')
        feqp1.write('%8d'%(ib_wfn + 1) + '%15.9f'%edft + '%15.9f'%eqp1 + '\n')

        fsig_scissors.write('%4d'%(ib_wfn + 1) + '%15.9f'%edft + '%15.9f'%edft + '            0.0' + '            0.0' + '            0.0' + '            0.0' + '            0.0' +  '%15.9f'%eqp0 + '%15.9f'%eqp1 + '            1.0' + '\n' )
    feqp0.close()
    feqp1.close()
    fsig_scissors.close()

  return

#endregion

#region classes
#endregion