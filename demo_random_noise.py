#!/usr/bin/env python3
import numpy as np
import util
import matplotlib.pyplot as plt

np.random.seed(1)

nsample = 2000
kmax = 99
nkx = 2*kmax+1
nky = kmax+1
n = 2*(kmax+1)
ii, jj = np.mgrid[0:n, 0:n]

# power_law = 1
# Pk = lambda k: k**((power_law-1)/2)
Pk = lambda k: np.exp(-(k-7)**2/10) + 0.3*np.exp(-((k-30)**2/20))
amplitude = 5

plt.switch_backend('Agg')
fig, ax = plt.subplots(1, 3, figsize=(12, 4))

#generating noise
noise = np.zeros((n, n, nsample))
for m in range(nsample):
  noise[:, :, m] = util.gaussian_random_field(Pk, n)
  noise[:, :, m] = (noise[:, :, m]-np.mean(noise[:, :, m]))*amplitude/np.std(noise[:, :, m])

#correlation map
ic, jc = int(n/2), int(n/2)
dist = np.sqrt((ii-ic)**2+(jj-jc)**2)
noise_mean = np.mean(noise, axis=2)
pert = noise.copy()
for m in range(nsample):
  pert[:, :, m] -= noise_mean
pert0 = pert[ic, jc, :]
var0 = np.sum(pert0**2) / (nsample-1)
var = np.sum(pert**2, axis=2) / (nsample-1)
cov = np.sum(pert * np.tile(pert0, (n, n, 1)), axis=2) / (nsample-1)
corr = cov / np.sqrt(var*var0)

# in terms of distance
dist1 = np.arange(0, n/2, 0.5)
corr1 = np.zeros(dist1.size)
for i in range(dist1.size):
  corr1[i] = np.mean(corr[np.where(abs(dist-dist1[i])<=0.5)])

##power spectrum
wn, pwr = util.pwrspec2d(noise)

ax[0].contourf(ii, jj, noise[:, :, 0], np.arange(-25, 25, 1), cmap='seismic')
ax[0].set_aspect('equal', 'box')
ax[0].set_xlabel('x', fontsize=15)
ax[0].set_ylabel('y', fontsize=15)
ax[0].set_title('error field', fontsize=20)

# ax[1].contourf(ii, jj, corr, np.arange(-1, 1, 0.1))
ax[1].plot(dist1, corr1, color='black')
ax[1].set_xlim(-1, n/2)
ax[1].set_ylim(-1, 1)
ax[1].set_xlabel('distance', fontsize=15)
ax[1].set_title('correlation', fontsize=20)

ax[2].loglog(wn, np.mean(pwr, axis=1), color='black')
ax[2].set_xlabel('wavenumber', fontsize=15)
ax[2].set_ylim(1e-4, 1e2)
ax[2].set_title('spectrum', fontsize=20)

fig.tight_layout()
plt.savefig('1.pdf')
