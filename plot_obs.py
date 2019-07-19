#!/usr/bin/env python3
import numpy as np
import util
import matplotlib.pyplot as plt

workdir = "/glade/scratch/mying/qgmodel_enkf/multiscale/"         #working directory
kmax = 63
nz = 2
t = 1

nx = 2*(kmax+1)
ny = 2*(kmax+1)
lv = 0  #vertical level of obs
obs_thin = 3

plt.switch_backend('Agg')
fig = plt.figure(figsize=(10, 5))

psik = util.read_field(workdir+"truth/{:05d}.bin".format(t), 2*kmax+1, kmax+1, nz)
xt = util.spec2grid(util.psi2temp(psik))[:, :, lv]
xt = np.roll(np.roll(xt, -40, axis=0), 60, axis=1)
ax = plt.subplot(121)
c = ax.contourf(xt.T, np.arange(-40, 40, 2), cmap='seismic')
plt.colorbar(c)
ax.set_aspect('equal', 'box')
ax.set_xlim(0, nx)
ax.set_ylim(0, ny)
ax.set_xticks(np.arange(0, nx+1, 20))
ax.set_yticks(np.arange(0, ny+1, 20))
ax.tick_params(labelsize=10)

xobs = xt + np.random.normal(0, 1, (nx, ny))
ax = plt.subplot(122)
obsmax = 40
obsmin = -40
cmap = [plt.cm.seismic(x) for x in np.linspace(0, 1, round(obsmax-obsmin)+1)]
for i in range(0, nx, obs_thin):
  for j in range(0, ny, obs_thin):
    ind = max(min(int(round(xobs[i, j] - obsmin)), int(round(obsmax-obsmin))), 0)
    ax.scatter(i, j, s=10, c=[cmap[ind][0:3]])
ax.set_aspect('equal', 'box')
ax.set_xlim(0, nx)
ax.set_ylim(0, ny)
ax.set_xticks(np.arange(0, nx+1, 20))
ax.set_yticks(np.arange(0, ny+1, 20))
ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig('1.pdf')
