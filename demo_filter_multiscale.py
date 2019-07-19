#!/usr/bin/env python3
##filter program: demo multiscale algorithm
import numpy as np
import util
import data_assimilation as da
import sys
import matplotlib.pyplot as plt

if len(sys.argv) != 8:
  print("usage: {} workdir casename kmax nz nens t run_displacement".format(sys.argv[0]))
  exit()

workdir = sys.argv[1]         #working directory
casename = sys.argv[2]        #casename for ensemble filtering experiment
kmax = int(sys.argv[3])       #maximum wavenumber (size of model n = 2*(kmax+1))
nz = int(sys.argv[4])         #number of vertical levels
nens = int(sys.argv[5])       #ensemble size
t = int(sys.argv[6])          #time step
run_displacement = int(sys.argv[7])

convertor = util.psi2temp
obs_ind = 7 #0=x, 1=y, 2=z, 3=u, 4=v, 5=psi, 6=zeta, 7=temp
obs_err = 1
krange = (3, 8, 20)
localize_cutoff = (24, 16, 10)
# krange = (1,)
# localize_cutoff = (16,)
clevel = (8, 5, 3)

nkx = 2*kmax+1
nky = kmax+1
nx = 2*(kmax+1)
ny = 2*(kmax+1)

###read truth
filename = workdir+'/truth/{:05d}.bin'.format(t)
psik = util.read_field(filename, nkx, nky, nz)
xt = util.spec2grid(convertor(psik))

###read prior ensemble
x = np.zeros((nens, nx, ny, nz))
for m in range(nens):
  filename = workdir+'/'+casename+'/'+'{:04d}'.format(m+1)+'/f_{:05d}.bin'.format(t)
  psik = util.read_field(filename, nkx, nky, nz)
  x[m, :, :, :] = util.spec2grid(convertor(psik))
xb = x.copy()

# util.output_ens('1.nc', xb[:, :, :, :])

###read obs
obsfile = workdir+'/obs/{:05d}'.format(t)
dat = np.loadtxt(obsfile)
nobs, nf = dat.shape
obs_loc = dat[:, 0:3] - 1.0
obs = dat[:, obs_ind]

##separate state variable into scale bands
ns = len(krange)
xs = np.zeros((nens, ns, nx, ny, nz))
xts = np.zeros((ns, nx, ny, nz))  ##truth
for s in range(ns):
  xts[s, :, :, :] = util.spec2grid(util.spec_bandpass(util.grid2spec(xt), krange, s))
  for m in range(nens):
    xs[m, s, :, :, :] = util.spec2grid(util.spec_bandpass(util.grid2spec(x[m, :, :, :]), krange, s))
xbs = xs.copy()  ##original prior copy
xs1 = xs.copy()
xs2 = xs.copy()

##run filter
u = np.zeros((nens, ns, nx, ny))
v = np.zeros((nens, ns, nx, ny))
for s in range(ns):
  print('running EnSRF for scale {}'.format(s))
  xs1[:, s, :, :, :] = xs[:, s, :, :, :]  #save a copy of the prior during iteration
  xs = da.EnSRF(xs, obs_loc, obs, obs_err, localize_cutoff[s], s)  #run DA step of iteration
  if run_displacement == 1:  #run displacement step of iteration
    if s < ns-1:
      print('aligning members')
      for m in range(nens):
        # print('aligning member {:04d}'.format(m+1))
        u[m, s, :, :], v[m, s, :, :] = da.optical_flow_HS(xs1[m, s, :, :, 0], xs[m, s, :, :, 0], nlevel=5)
        for z in range(nz):
          for i in range(s+1, ns):
            xs[m, i, :, :, z] = util.warp(xs[m, i, :, :, z], -u[m, s, :, :], -v[m, s, :, :])
        xs2[:, s+1, :, :, :] = xs[:, s+1, :, :, :] #save a copy of the aligned prior
print('filtering complete')
xas = xs.copy() #final analysis

# util.output_ens('2.nc', xa[:, :, :, :])

def set_axis(ax, title):
  ax.set_aspect('equal', 'box')
  ax.set_xlim(0, nx)
  ax.set_ylim(0, ny)
  ax.set_xticks(np.arange(0, nx+1, 50))
  ax.set_yticks(np.arange(0, ny+1, 50))
  ax.tick_params(labelsize=10)
  ax.set_title(title, fontsize=15)

#shift domain position for better plotting
xt = np.roll(np.roll(xt, -40, axis=0), 60, axis=1)
xb = np.roll(np.roll(xb, -40, axis=1), 60, axis=2)
xts = np.roll(np.roll(xts, -40, axis=1), 60, axis=2)
xbs = np.roll(np.roll(xbs, -40, axis=2), 60, axis=3)
xas = np.roll(np.roll(xas, -40, axis=2), 60, axis=3)
xs2 = np.roll(np.roll(xs2, -40, axis=2), 60, axis=3)
u = np.roll(np.roll(u, -40, axis=2), 60, axis=3)
v = np.roll(np.roll(v, -40, axis=2), 60, axis=3)

for m in range(0, 1):
  plt.switch_backend('Agg')
  fig, ax = plt.subplots(ns+1, 4, figsize=(4*4, 4*(ns+1)))
  ii, jj = np.mgrid[0:nx, 0:ny]
  lv = 0
  cmap = [plt.cm.jet(m) for m in np.linspace(0, 1, nens)]
  for s in range(ns):
    print('diagnostics for scale {}, ROI={}'.format(s, localize_cutoff[s]))
    ax[s, 0].contourf(ii, jj, xts[s, :, :, lv], np.arange(-40, 40, 2), cmap='seismic')
    ax[s, 0].contour(ii, jj, xts[s, :, :, lv], (-clevel[s], clevel[s]), linestyles='-', colors='k')
    set_axis(ax[s, 0], 'truth')

    ax[s, 1].contourf(ii, jj, np.mean(xbs[m:m+1, s, :, :, lv], axis=0), np.arange(-40, 40, 2), cmap='seismic')
    ax[s, 1].contour(ii, jj, np.mean(xbs[m:m+1, s, :, :, lv], axis=0), (-clevel[s], clevel[s]), linestyles='-', colors='r')
    ax[s, 1].contour(ii, jj, xts[s, :, :, lv], (-clevel[s], clevel[s]), linestyles='-', colors='k')
    set_axis(ax[s, 1], 'prior')
    print('  prior     rmse={:5.2f}, {:5.2f}'.format(util.rmse(xbs[m, s, :, :, lv], xts[s, :, :, lv]), util.rmse(np.mean(xbs[:, s, :, :, lv], axis=0), xts[s, :, :, lv])))

    ax[s, 2].contourf(ii, jj, np.mean(xs2[m:m+1, s, :, :, lv], axis=0), np.arange(-40, 40, 2), cmap='seismic')
    ax[s, 2].contour(ii, jj, np.mean(xs2[m:m+1, s, :, :, lv], axis=0), (-clevel[s], clevel[s]), linestyles='-', colors='r')
    ax[s, 2].contour(ii, jj, xts[s, :, :, lv], (-clevel[s], clevel[s]), linestyles='-', colors='k')
    qv = ax[s, 2].quiver(ii[::3, ::3], jj[::3, ::3], u[m, s, ::3, ::3], v[m, s, ::3, ::3], scale=100, headwidth=4)
    set_axis(ax[s, 2], 'aligned prior')

    ax[s, 3].contourf(ii, jj, np.mean(xas[m:m+1, s, :, :, lv], axis=0), np.arange(-40, 40, 2), cmap='seismic')
    ax[s, 3].contour(ii, jj, np.mean(xas[m:m+1, s, :, :, lv], axis=0), (-clevel[s], clevel[s]), linestyles='-', colors='r')
    ax[s, 3].contour(ii, jj, xts[s, :, :, lv], (-clevel[s], clevel[s]), linestyles='-', colors='k')
    set_axis(ax[s, 3], 'posterior')
    print('  posterior rmse={:5.2f}, {:5.2f}'.format(util.rmse(xas[m, s, :, :, lv], xts[s, :, :, lv]), util.rmse(np.mean(xas[:, s, :, :, lv], axis=0), xts[s, :, :, lv])))

  xa = np.sum(xas, axis=1)
  print('full prior     state rmse = {:5.2f}, sprd = {:5.2f}'.format(util.rmse(np.mean(xb[:, :, :, lv], axis=0), xt[:, :, lv]),util.sprd(xb[:, :, :, lv])))
  print('full posterior state rmse = {:5.2f}, sprd = {:5.2f}'.format(util.rmse(np.mean(xa[:, :, :, lv], axis=0), xt[:, :, lv]),util.sprd(xa[:, :, :, lv])))

  # plt.savefig('{:02d}.pdf'.format(m))
  plt.savefig('1.pdf')
  plt.close()
