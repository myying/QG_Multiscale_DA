#!/usr/bin/env python3
###read in ensemble state from bin files and plot error spectrum
import numpy as np
import util
import sys

if len(sys.argv) != 5:
  print("usage: {} workdir kmax nz nt".format(sys.argv[0]))
  exit()

workdir = sys.argv[1]
kmax = int(sys.argv[2])
nz = int(sys.argv[3])
nt = int(sys.argv[4])

nkx = 2*kmax+1
nky = kmax+1
nx = 2*(kmax+1)
ny = 2*(kmax+1)
convert_var = util.psi2temp

nup = int(np.ceil((max(nx, ny)+1)/2))
pwr = np.zeros((nup, nz, nt))

for t in range(nt):
  filename = workdir+'/truth/{:05d}.bin'.format(t+1)
  psik = util.read_field(filename, nkx, nky, nz)
  var_truth = util.spec2grid(convert_var(psik))
  wn, pwr[:, :, t] = util.pwrspec2d(var_truth)  #reference spectrum from truth

##save data
np.save(workdir+'/truth/spectrum', pwr)
