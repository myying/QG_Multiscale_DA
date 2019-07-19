#!/usr/bin/env python3
import numpy as np
import util
import sys

if len(sys.argv) != 8:
  print("usage: {} workdir casename kmax nz nens amplitude power_law".format(sys.argv[0]))
  exit()

workdir = sys.argv[1]
casename = sys.argv[2]
kmax = int(sys.argv[3])
nz = int(sys.argv[4])
nens = int(sys.argv[5])
amplitude = float(sys.argv[6])
power_law = float(sys.argv[7])

n = 2*(kmax+1)

#read truth initial condition
psik = util.read_field(workdir+'/initial_condition.bin', 2*kmax+1, kmax+1, nz)
psi = util.spec2grid(psik)

Pk = lambda k: k**((power_law-1)/2)

for m in range(nens):
  psi1 = psi.copy()
  noise = util.gaussian_random_field(Pk, n)
  noise = noise*amplitude/np.std(noise)
  for z in range(nz):
    psi1[:, :, z] += noise
  psik1 = util.grid2spec(psi1)
  util.write_field(workdir+'/'+casename+'/{:04d}/f_00001.bin'.format(m+1), psik1)

