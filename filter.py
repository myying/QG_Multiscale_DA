#!/usr/bin/env python3
##filter program
import numpy as np
import util
import data_assimilation as da
import sys

if len(sys.argv) != 9:
  print("usage: {} workdir casename kmax nz nens t run_displacement run_multiscale".format(sys.argv[0]))
  exit()

workdir = sys.argv[1]         #working directory
casename = sys.argv[2]        #casename for ensemble filtering experiment
kmax = int(sys.argv[3])       #maximum wavenumber (size of model n = 2*(kmax+1))
nz = int(sys.argv[4])         #number of vertical levels
nens = int(sys.argv[5])       #ensemble size
t = int(sys.argv[6])          #time step
run_displacement = int(sys.argv[7])
run_multiscale = int(sys.argv[8])

print('analysis cycle #{:05d}'.format(t))

convertor = util.psi2temp
convertor_back = util.temp2psi
obs_ind = 7 #0=x, 1=y, 2=z, 3=u, 4=v, 5=psi, 6=zeta, 7=temp
obs_err = 1
if run_multiscale == 1:
  krange = (3, 8, 20)
  localize_cutoff = (24, 16, 10)
  if nens == 5:
    localize_cutoff = (12, 8, 5)
  if nens == 10:
    localize_cutoff = (18, 12, 7)
  if nens == 20:
    localize_cutoff = (24, 16, 10)
  if nens > 20:
    localize_cutoff = (30, 22, 15)
else:
  krange = (1,)
  localize_cutoff = (16,)
  if nens == 5:
    localize_cutoff = (8,)
  if nens == 10:
    localize_cutoff = (12,)
  if nens == 20:
    localize_cutoff = (16,)
  if nens > 20:
    localize_cutoff = (22,)
print(localize_cutoff)

nkx = 2*kmax+1
nky = kmax+1
nx = 2*(kmax+1)
ny = 2*(kmax+1)

###read prior ensemble
x = np.zeros((nens, nx, ny, nz))
for m in range(nens):
  filename = workdir+'/'+casename+'/'+'{:04d}'.format(m+1)+'/f_{:05d}.bin'.format(t)
  psik = util.read_field(filename, nkx, nky, nz)
  x[m, :, :, :] = util.spec2grid(convertor(psik))
xb = x.copy()

###read obs
obsfile = workdir+'/obs/{:05d}'.format(t)
dat = np.loadtxt(obsfile)
nobs, nf = dat.shape
obs_loc = dat[:, 0:3] - 1.0
obs = dat[:, obs_ind]

##separate state variable into scale bands
ns = len(krange)
xs = np.zeros((nens, ns, nx, ny, nz))
for s in range(ns):
  for m in range(nens):
    xs[m, s, :, :, :] = util.spec2grid(util.spec_bandpass(util.grid2spec(x[m, :, :, :]), krange, s))

##run filter
u = np.zeros((nens, ns, nx, ny))
v = np.zeros((nens, ns, nx, ny))
for s in range(ns):
  xs1 = xs.copy()
  # print('running EnSRF for scale {}'.format(s))
  xs = da.EnSRF(xs, obs_loc, obs, obs_err, localize_cutoff[s], s)
  if run_displacement == 1:
    if s < ns-1:
      for m in range(nens):
        # print('aligning member {:04d}'.format(m+1))
        u[m, s, :, :], v[m, s, :, :] = da.optical_flow_HS(xs1[m, s, :, :, 0], xs[m, s, :, :, 0], nlevel=5)
        for z in range(nz):
          for i in range(s+1, ns):
            xs[m, i, :, :, z] = util.warp(xs[m, i, :, :, z], -u[m, s, :, :], -v[m, s, :, :])
xas = xs.copy()

##sum scales back to full state
xa = np.sum(xas, axis=1)

###adaptive inflation
hxb = np.zeros((nens, nobs))
hxa = np.zeros((nens, nobs))
for m in range(nens):
  for n in range(nobs):
    hxb[m, n] = util.interp3d(xb[m, :, :, :], obs_loc[n, :])
    hxa[m, n] = util.interp3d(xa[m, :, :, :], obs_loc[n, :])
hxbm = np.mean(hxb, axis=0)
hxam = np.mean(hxa, axis=0)
amb = hxam - hxbm
oma = obs - hxam
vara = np.sum((hxa - np.tile(hxam, (nens, 1)))**2, axis=0)/(nens-1)
infl = np.sqrt(max(1.0, np.sum(amb*oma)/np.sum(vara)))
infl = min(2.0, infl)
print('adaptive inflation lambda={:5.3f}'.format(infl))

xam = np.mean(xa, axis=0)
for m in range(nens):
  xa[m, :, :, :] = xam + infl*(xa[m, :, :, :] - xam)

##output aligned prior ensemble
xb_align = xb.copy()
for m in range(nens):
  for z in range(nz):
    xb_align[m, :, :, z] = util.warp(xb[m, :, :, z], -np.sum(u[m, :, :, :], axis=0), -np.sum(v[m, :, :, :], axis=0))

##output posterior ensemble
for m in range(nens):
  ##preserve k=0 in psi
  psik0 = util.read_field(workdir+'/'+casename+'/'+'{:04d}'.format(m+1)+'/f_{:05d}.bin'.format(t), nkx, nky, nz)
  psik = convertor_back(util.grid2spec(xa[m, :, :, :]))
  psik[kmax, 0, :] = psik0[kmax, 0, :]
  filename = workdir+'/'+casename+'/'+'{:04d}'.format(m+1)+'/{:05d}.bin'.format(t)
  util.write_field(filename, psik)

  psik = convertor_back(util.grid2spec(xb_align[m, :, :, :]))
  psik[kmax, 0, :] = psik0[kmax, 0, :]
  filename = workdir+'/'+casename+'/'+'{:04d}'.format(m+1)+'/fa_{:05d}.bin'.format(t)
  util.write_field(filename, psik)
