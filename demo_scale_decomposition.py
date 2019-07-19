#!/usr/bin/env python3
import numpy as np
import util
import matplotlib.pyplot as plt

kmax = 63
nz = 2
lv = 0
psik = util.read_field('/glade/scratch/mying/qgmodel_enkf/multiscale/truth/00001.bin', 2*kmax+1, kmax+1, nz)
var = util.spec2grid(util.psi2temp(psik))
wn, pwr = util.pwrspec2d(var)
specsmth = 1.0
krange = (3, 8, 20)
ns = len(krange)
n = 2*(kmax+1)
ii, jj = np.mgrid[0:n, 0:n]

plt.switch_backend('Agg')
fig, ax = plt.subplots(2, ns+1, figsize=(4*(ns+1), 8))

pwr_s = np.zeros((ns, wn.size, nz))
var_s = np.zeros((ns, n, n, nz))
for s in range(ns):
  var_s[s, :, :, :] = util.spec2grid(util.spec_bandpass(util.grid2spec(var), krange, s))
  wn, pwr_s[s, :, :] = util.pwrspec2d(var_s[s, :, :, :])

print(np.sum(np.sum(var_s, axis=0) - var)) ###check if sum back to var

for s in range(ns):
  for i in range(ns):
    ax[0, s].loglog(wn, util.smooth_spec(wn, pwr_s[i, :, lv], specsmth), color=[.7, .7, .7])
  ax[0, s].loglog(wn, util.smooth_spec(wn, pwr_s[s, :, lv], specsmth), color='k', linewidth=2)
  ax[0, s].set_xlim(1, 100)
  ax[0, s].set_ylim(1e-3, 1e2)

  ax[1, s].contourf(ii, jj, var_s[s, :, :, lv], np.arange(-40, 40, 2), cmap='seismic')
  ax[1, s].set_aspect('equal', 'box')


plt.savefig('1.pdf')
