###prototype data assimilation algorithms
import numpy as np
import util

###Ensemble square root filter, a serial variant of deterministic EnKF (Whitaker and Hamill 2002)
def EnSRF(x, obs_loc, obs, obserr, localize_cutoff, s):
  nens, ns, ni, nj, nk = x.shape
  nobs = obs.size
  ii, jj, kk = np.mgrid[0:ni, 0:nj, 0:nk]

  for n in range(nobs):
    xm = np.mean(x[:, s, :, :, :], axis=0)
    xp = x[:, s, :, :, :] - np.tile(xm, (nens, 1, 1, 1))
    yo = obs[n]
    io = obs_loc[n, 0]
    jo = obs_loc[n, 1]
    ko = obs_loc[n, 2]
    hx = np.zeros(nens)
    for m in range(nens):
      hx[m] = util.interp3d(np.sum(x[m, :, :, :, :], axis=0), obs_loc[n, :])
    hxm = np.mean(hx)
    hxp = hx - hxm

    innov = yo - hxm
    varo = obserr**2
    varb = np.sum(hxp**2) / (nens - 1)
    cov = np.sum(xp * np.tile(hxp, (nk, nj, ni, 1)).T, axis=0) / (nens - 1)

    dist = np.sqrt(np.minimum(np.abs(ii-io),ni-np.abs(ii-io))**2+np.minimum(np.abs(jj-jo),nj-np.abs(jj-jo))**2)

    loc = local_GC(dist, localize_cutoff)
    gain = loc * cov / (varo + varb)
    srf = 1.0 / (1.0 + np.sqrt(varo / (varo + varb)))

    xm = xm + gain * innov
    for m in range(nens):
      xp[m, :, :, :] = xp[m, :, :, :] - srf * gain * hxp[m]
    x[:, s, :, :, :] = xp + np.tile(xm, (nens, 1, 1, 1))

    # print('{:-6d} ({:6.1f},{:6.1f},{:4.1f}) yo={:7.2f} hxm={:7.2f} varb={:9.2f}'.format(n+1, obs_loc[n, 0], obs_loc[n, 1], obs_loc[n, 2], obs[n], hxm, varb))

  return x


####Adaptive covariance relaxation method (Ying and Zhang 2015, QJRMS)
###ensures ensemble spread is large enough during cycling
def adaptive_relax_coef(xb, xa, obs_loc, obs, obserr):
  nobs = obs.size
  nens, nx, ny, nz = xb.shape
  hxb = np.zeros((nens, nobs))
  hxa = np.zeros((nens, nobs))
  for n in range(nobs):
    for m in range(nens):
      hxb[m, n] = util.interp3d(xb[m, :, :, :], obs_loc[n, :])
      hxa[m, n] = util.interp3d(xa[m, :, :, :], obs_loc[n, :])
  hxbm = np.mean(hxb, axis=0)
  hxam = np.mean(hxa, axis=0)
  omaamb = np.sum((obs-hxam)*(hxam-hxbm))/nobs
  omb2 = np.sum((obs-hxbm)**2)/nobs
  amb2 = np.sum((hxam-hxbm)**2)/nobs
  varo = obserr**2
  varb = np.sum(np.std(hxb, axis=0)**2)/nobs
  vara = np.sum(np.std(hxa, axis=0)**2)/nobs
  beta = np.sqrt(varb/vara)
  lamb = np.sqrt(max(0.0, (omb2-varo-amb2)/vara))
  alpha = (lamb-1)/(beta-1)
  if alpha > 2:
    alpha = 2
  if alpha < -1:
    alpha = -1
  if beta < 1:
    alpha = 0
  print('adaptive covariance relaxation:')
  print('var_b = {:7.2f}, var_a = {:7.2f}'.format(varb, vara))
  print('omb2  = {:7.2f}, amb2  = {:7.2f}'.format(omb2, amb2))
  print('beta  = {:7.5f}, lambda = {:7.5f}'.format(beta, lamb))
  print('alpha = {:7.5f}'.format(alpha))
  return alpha

####alignment techniques
###compute optical flow in observed space using hierarchical HS algorithm
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


###localization functions
###Gaspari Cohn 1999:
def local_GC(dist, cutoff):
  loc = np.zeros(dist.shape)
  if cutoff>0:
    r = dist / (cutoff / 2)
    loc1 = (((-0.25*r + 0.5)*r + 0.625)*r - 5.0/3.0) * r**2 + 1
    ind1 = np.where(dist<cutoff/2)
    loc[ind1] = loc1[ind1]
    r[np.where(r==0)] = 1e-10
    loc2 = ((((r/12.0 - 0.5)*r + 0.625)*r + 5.0/3.0)*r - 5.0)*r + 4 - 2.0/(3.0*r)
    ind2 = np.where(np.logical_and(dist>=cutoff/2, dist<cutoff))
    loc[ind2] = loc2[ind2]
  else:
    loc = np.ones(dist.shape)
  return loc

###Gaussian localization function:
def local_Gauss(dist, cutoff):
  loc = np.exp(-dist**2 / (cutoff / 3.5)**2 / 2)
  return loc
