#!/usr/bin/env python3
import numpy as np
import util
import sys

if len(sys.argv) != 5:
  print("usage: {} kmax nz amplitude power_law".format(sys.argv[0]))
  exit()

kmax = int(sys.argv[1])
nz = int(sys.argv[2])
amplitude = float(sys.argv[3])
power_law = float(sys.argv[4])

n = 2*(kmax+1)
Pk = lambda k: k**((power_law-1)/2)
noise = util.gaussian_random_field(Pk, n)
print(np.std(noise))
noise = noise*amplitude/np.std(noise)
ic = np.zeros((n, n, nz))
for z in range(nz):
  ic[:, :, z] = noise

ick = util.grid2spec(ic)

util.write_field('initial_condition.bin', ick)

