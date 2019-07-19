#!/usr/bin/env python3
###read in ensemble state from bin files and plot error spectrum
import numpy as np
import util
import sys

if len(sys.argv) != 8:
  print("usage: {} workdir casename kmax nz nens nt ens_type".format(sys.argv[0]))
  exit()

workdir = sys.argv[1]
casename = sys.argv[2]
kmax = int(sys.argv[3])
nz = int(sys.argv[4])
nens = int(sys.argv[5])
nt = int(sys.argv[6])
ens_type = int(sys.argv[7])  #1: prior ensemble, 2: posterior ensemble

nkx = 2*kmax+1
nky = kmax+1
nx = 2*(kmax+1)
ny = 2*(kmax+1)
convert_var = util.psi2temp

nup = int(np.ceil((max(nx, ny)+1)/2))
err = np.zeros((nup, nz, nt))
pwr_ens = np.zeros((nup, nz, nens, nt))

for t in range(0, nt):
  filename = workdir+'/truth/{:05d}.bin'.format(t+1)
  psik = util.read_field(filename, nkx, nky, nz)
  var_truth = util.spec2grid(convert_var(psik))

  var_ens = np.zeros((nens, nx, ny, nz))
  if ens_type == 1:
    name = 'f_{:05d}'.format(t+1)
    disp = 'prior'
  else:
    name = '{:05d}'.format(t+1)
    disp = 'post'
  for m in range(nens):
    filename = workdir+'/'+casename+'/'+'{:04d}'.format(m+1)+'/'+name+'.bin'
    psik = util.read_field(filename, nkx, nky, nz)
    var_ens[m, :, :, :] = util.spec2grid(convert_var(psik)) #ensemble members
  var_mean = np.mean(var_ens, axis=0) #ensemble mean

  for m in range(nens):
    wn, pwr_ens[:, :, m, t] = util.pwrspec2d(var_ens[m, :, :, :] - var_mean)

  wn, err[:, :, t] = util.pwrspec2d(var_mean - var_truth) #ensemble mean error spectrum

##save data
np.save(workdir+'/'+casename+'/spectrum_'+disp+'_err', err)
np.save(workdir+'/'+casename+'/spectrum_'+disp+'_ens', pwr_ens)
