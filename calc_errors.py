#!/usr/bin/env python3
###read in ensemble state from bin files and calculate errors/spread
import numpy as np
import util
import sys

if len(sys.argv) != 8:
  print("usage: {} workdir casename kmax nz nens cp ncycle".format(sys.argv[0]))
  exit()

workdir = sys.argv[1]
casename = sys.argv[2]
kmax = int(sys.argv[3])
nz = int(sys.argv[4])
nens = int(sys.argv[5])
cp = int(sys.argv[6])
ncycle = int(sys.argv[7])

nkx = 2*kmax+1
nky = kmax+1
nx = 2*(kmax+1)
ny = 2*(kmax+1)
convert_var = util.psi2temp

rmse = np.zeros((ncycle, 2))
sprd = np.zeros((ncycle, 2))

for n in range(ncycle):
  name = '{:05d}.bin'.format(n*cp+1)

  filename = workdir+'/truth/'+name
  psik = util.read_field(filename, nkx, nky, nz)
  var_truth = util.spec2grid(convert_var(psik))

  var_ens = np.zeros((nens, nx, ny, nz, 2))
  for m in range(nens):
    filename = workdir+'/'+casename+'/{:04d}'.format(m+1)+'/f_'+name
    psik = util.read_field(filename, nkx, nky, nz)
    var = util.spec2grid(convert_var(psik)) #ensemble members
    var_ens[m, :, :, :, 0] = var
    filename = workdir+'/'+casename+'/{:04d}'.format(m+1)+'/'+name
    psik = util.read_field(filename, nkx, nky, nz)
    var = util.spec2grid(convert_var(psik)) #ensemble members
    var_ens[m, :, :, :, 1] = var
  var_mean = np.mean(var_ens, axis=0) #ensemble mean

  for i in range(2):
    rmse[n, i] = np.sqrt(np.mean((var_mean[:, :, :, i] - var_truth)**2))
    epvar = np.sum((var_ens[:, :, :, :, i] - np.tile(var_mean[:, :, :, i], (nens, 1, 1, 1)))**2, axis=0)/(nens-1)
    sprd[n, i] = np.sqrt(np.mean(epvar))

##save data
np.save(workdir+'/'+casename+'/rmse', rmse)
np.save(workdir+'/'+casename+'/sprd', sprd)
