#!/usr/bin/env python3
import numpy as np
import util
import sys

if len(sys.argv) != 6:
  print("usage: {} workdir kmax nz nt obs_thin".format(sys.argv[0]))
  exit()

workdir = sys.argv[1]
kmax = int(sys.argv[2])
nz = int(sys.argv[3])
nt = int(sys.argv[4])
obs_thin = int(sys.argv[5])
obs_z = (0,)

n = 2*(kmax+1)
power_law = 1
Pk = lambda k: k**((power_law-1)/2)
psi_err = 0.5
uv_err = 1
temp_err = 1
zeta_err = 30
smth = 0

nx = 2*(kmax+1)
ny = 2*(kmax+1)
xx, yy, zz = np.mgrid[1:nx+1, 1:ny+1, 1:nz+1]

for t in range(nt):
  psik = util.read_field(workdir+'/truth/{:05d}.bin'.format(t+1), 2*kmax+1, kmax+1, nz)
  psi0 = util.spec2grid(psik)
  u0 = util.spec2grid(util.psi2u(psik))
  v0 = util.spec2grid(util.psi2v(psik))
  temp0 = util.spec2grid(util.psi2temp(psik))
  zeta0 = util.spec2grid(util.psi2zeta(psik))

  #add noise
  noise = util.gaussian_random_field(Pk, n)
  psi = psi0 + psi_err*np.tile(noise, (nz, 1, 1)).T
  noise = util.gaussian_random_field(Pk, n)
  u = u0 + uv_err*np.tile(noise, (nz, 1, 1)).T
  noise = util.gaussian_random_field(Pk, n)
  v = v0 + uv_err*np.tile(noise, (nz, 1, 1)).T
  noise = util.gaussian_random_field(Pk, n)
  temp = temp0 + temp_err*np.tile(noise, (nz, 1, 1)).T
  noise = util.gaussian_random_field(Pk, n)
  zeta = zeta0 + zeta_err*np.tile(noise, (nz, 1, 1)).T

  #thinning
  xx1 = xx[::obs_thin, ::obs_thin, obs_z]
  yy1 = yy[::obs_thin, ::obs_thin, obs_z]
  zz1 = zz[::obs_thin, ::obs_thin, obs_z]
  psi1 = psi[::obs_thin, ::obs_thin, obs_z]
  u1 = u[::obs_thin, ::obs_thin, obs_z]
  v1 = v[::obs_thin, ::obs_thin, obs_z]
  temp1 = temp[::obs_thin, ::obs_thin, obs_z]
  zeta1 = zeta[::obs_thin, ::obs_thin, obs_z]
  nx1, ny1, nz1 = xx1.shape

  #smoothing
  psi1 = util.smooth(psi1, smth)
  u1 = util.smooth(u1, smth)
  v1 = util.smooth(v1, smth)
  temp1 = util.smooth(temp1, smth)
  zeta1 = util.smooth(zeta1, smth)

  f = open(workdir+'/obs/{:05d}'.format(t+1), 'w')
  for z in range(nz1):
    for y in range(ny1):
      for x in range(nx1):
        f.write('{:7.2f} {:7.2f} {:5.2f} {:12.5f} {:12.5f} {:12.5f} {:12.5f} {:12.5f}\n'.format(xx1[x,y,z], yy1[x,y,z], zz1[x,y,z], u1[x,y,z], v1[x,y,z], psi1[x,y,z], zeta1[x,y,z], temp1[x,y,z]))

  f.close()
