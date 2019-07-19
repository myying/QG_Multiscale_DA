#!/usr/bin/env python3
###read in ensemble state from bin files and plot error spectrum
import numpy as np
import util
import matplotlib.pyplot as plt

np.random.seed(1)
plt.switch_backend('Agg')
ax = plt.subplot(111)
specsmth = 1.1

kmax = 255
nkx = 2*kmax+1
nky = kmax+1
n = 2*(kmax+1)

###truth
power_law0 = -2
amplitude0 = 1e3
Pk0 = lambda k: k**((power_law0-1)/2)
truth = np.zeros((n, n, 1))
truth[:, :, 0] = util.gaussian_random_field(Pk0, n) * amplitude0
wn, pwr = util.pwrspec2d(truth)
ax.loglog(wn, util.smooth_spec(wn, pwr, specsmth), color='black', linewidth=2)

###noise
power_law = 1
amplitude = 3e2
Pk = lambda k: k**((power_law-1)/2)
# Pk = lambda k: np.exp(-(k-30)**2/100**2)

noise = np.zeros((n, n, 1))
noise[:, :, 0] = util.gaussian_random_field(Pk, n) * amplitude
obs = truth + noise

colors = ('red', 'blue', [.3, .7, .3], [.8, .8, .3], [.7, .3, .7])
###demonstrate thinning of obs
# intv = (1, 2, 4, 8)
# for i in range(len(intv)):
  # wn, pwr = util.pwrspec2d(noise[::intv[i], ::intv[i], :])
  # ax.loglog(wn, util.smooth_spec(wn, pwr, specsmth), color=colors[i], linewidth=2)

###demonstrate smoothing of obs
intv = 2
smth = (1, 2, 4, 8)
for i in range(len(smth)):
  obs_smth = util.smooth(obs[::intv, ::intv, :], smth[i])
  obs_error = obs_smth - truth[::intv, ::intv, :]
  wn, pwr = util.pwrspec2d(obs_error)
  ax.loglog(wn, util.smooth_spec(wn, pwr, specsmth), color=colors[i], linewidth=2)

###demonstrate downscaling effect
# ns = (n, 2*n, 4*n)
# for i in range(len(ns)):
#   noise_ds = np.zeros((ns[i], ns[i], 1))
#   noise_ds[:, :, 0] = util.regrid(noise[:, :, 0], ns[i], ns[i])
#   wn, pwr = util.pwrspec2d(noise_ds)
#   ax.loglog(wn, util.smooth_spec(wn, pwr, specsmth), color=colors[i], linewidth=2)

ax.set_xlim(1, 1000)
ax.set_ylim(1e0, 1e5)
ax.tick_params(labelsize=12)
ax.set_xlabel('wavenumber', fontsize=12)

plt.savefig('1.pdf')
