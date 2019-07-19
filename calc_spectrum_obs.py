#!/usr/bin/env python3
###read in ensemble state from bin files and plot error spectrum
import numpy as np
import util
import sys

if len(sys.argv) != 6:
  print("usage: {} workdir obsdir kmax nz nt".format(sys.argv[0]))
  exit()

workdir = sys.argv[1]
obsdir = sys.argv[2]
kmax = int(sys.argv[3])
nz = int(sys.argv[4])
nt = int(sys.argv[5])

convert_var = util.psi2temp
obs_ind = 7  # x, y, z, u, v, psi, zeta, temp
             # 0, 1, 2, 3, 4,   5,    6,    7
obs_thin = 3
obs_lv = 0

nx = 2*(kmax+1)
ny = 2*(kmax+1)

nxob = int(np.ceil(nx/obs_thin))
nyob = int(np.ceil(ny/obs_thin))

nupob = int(np.ceil((max(nxob, nyob)+1)/2))
obserr_pwr = np.zeros((nupob, 1, nt))

for t in range(nt):
  filename = workdir+'/truth/{:05d}.bin'.format(t+1)
  psik = util.read_field(filename, 2*kmax+1, kmax+1, nz)
  var_truth = util.spec2grid(convert_var(psik))
  obs_truth = var_truth[::obs_thin, ::obs_thin, obs_lv]

  dat = np.loadtxt(workdir+'/'+obsdir+'/{:05d}'.format(t+1))
  obs = np.reshape(dat[:, obs_ind], (nxob, nyob)).T

  obserr = np.zeros((nxob, nyob, 1))
  obserr[:, :, 0] = obs - obs_truth

  wn, obserr_pwr[:, :, t] = util.pwrspec2d(obserr)

##save data
np.save(workdir+'/'+obsdir+'/spectrum_obserr', obserr_pwr)
