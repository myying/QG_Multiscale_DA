#!/usr/bin/env python3
import numpy as np
import util
import matplotlib.pyplot as plt

plt.switch_backend('Agg')
plt.figure(figsize=(8, 8))

t = 1
kmax = 63
nz = 2
ni = 2*(kmax+1)
nj = 2*(kmax+1)
istart = 0
iend = ni-1
jstart = 0
jend = nj-1
smth = 0

def set_axis(ax, title):
  ax.set_aspect('equal', 'box')
  ax.set_xlim(istart, iend)
  ax.set_ylim(jstart, jend)
  ax.set_xticks(np.arange(istart, iend, 40))
  ax.set_yticks(np.arange(jstart, jend, 40))
  ax.set_title(title, fontsize=15)

workdir = '/glade/scratch/mying/qgmodel_enkf/multiscale'

psik = util.read_field(workdir+'/truth/{:05d}.bin'.format(t), kmax*2+1, kmax+1, nz)
xt = util.smooth(util.spec2grid(util.psi2temp(psik))[:, :, 0], smth)
xt = np.roll(np.roll(xt, -40, axis=0), 60, axis=1)

filename = workdir+'/noda/{:04d}/f_{:05d}.bin'.format(1, t)
psik = util.read_field(filename, kmax*2+1, kmax+1, nz)
xb = util.smooth(util.spec2grid(util.psi2temp(psik))[:, :, 0], smth)
xb = np.roll(np.roll(xb, -40, axis=0), 60, axis=1)

def optical_flow_HS(Im1, Im2, nlevel):
  ni, nj = Im1.shape
  u = np.zeros((ni, nj))
  v = np.zeros((ni, nj))
  for lev in range(nlevel, -1, -1):
    Im1warp = util.warp(Im1, -u, -v)
    Im1c = util.coarsen(Im1warp, lev)
    Im2c = util.coarsen(Im2, lev)

    niter = 100
    w1 = 100
    w2 = 0
    Ix = 0.5*(util.deriv_x(Im1c) + util.deriv_x(Im2c))
    Iy = 0.5*(util.deriv_y(Im1c) + util.deriv_y(Im2c))
    It = Im2c - Im1c
    du = np.zeros(Ix.shape)
    dv = np.zeros(Ix.shape)
    for k in range(niter):
      ubar2 = util.laplacian(du) + du
      vbar2 = util.laplacian(dv) + dv
      ubar1 = util.deriv_xx(du) + du
      vbar1 = util.deriv_yy(dv) + dv
      uxy = util.deriv_xy(du)
      vxy = util.deriv_xy(dv)
      du = (w1*ubar2 + w2*(ubar1+vxy))/(w1+w2) - Ix*((w1*(Ix*ubar2 + Iy*vbar2) + w2*((ubar1+vxy)*Ix + (vbar1+uxy)*Iy))/(w1+w2) + It)/(w1 + w2 + Ix**2 + Iy**2)
      dv = (w1*vbar2 + w2*(vbar1+uxy))/(w1+w2) - Iy*((w1*(Ix*ubar2 + Iy*vbar2) + w2*((ubar1+vxy)*Ix + (vbar1+uxy)*Iy))/(w1+w2) + It)/(w1 + w2 + Ix**2 + Iy**2)

    u += util.sharpen(du*2**lev, lev)
    v += util.sharpen(dv*2**lev, lev)
  return u, v

ii, jj = np.mgrid[0:ni, 0:nj]
ax = plt.subplot(2,2,1)
ax.contourf(ii, jj, xt, np.arange(-40, 40, 2), cmap='seismic')
ax.contour(ii, jj, xt, (-10, 10), colors='k', linestyles='solid')
set_axis(ax, '(a) target field')
ax = plt.subplot(2,2,2)
ax.contourf(ii, jj, xb, np.arange(-40, 40, 2), cmap='seismic')
# ax.contour(ii, jj, xb, (-10, 10), colors='gray', linestyles='solid')
ax.contour(ii, jj, xt, (-10, 10), colors='k', linestyles='solid')
set_axis(ax, '(b) source field')

ax = plt.subplot(2,2,3)
u, v = optical_flow_HS(xb, xt, 7)
ax.quiver(ii[::3, ::3], jj[::3, ::3], u[::3, ::3], v[::3, ::3], scale=100, headwidth=4)
set_axis(ax, '(c) displacement')

ax = plt.subplot(2,2,4)
xa = util.warp(xb, -u, -v)
ax.contourf(ii, jj, xa, np.arange(-40, 40, 2), cmap='seismic')
# ax.contour(ii, jj, xa, (-10, 10), colors='gray', linestyles='solid')
ax.contour(ii, jj, xt, (-10, 10), colors='k', linestyles='solid')
set_axis(ax, '(d) aligned field')


plt.savefig('1.pdf')
