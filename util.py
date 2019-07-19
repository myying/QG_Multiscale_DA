###Utility function for handling qgmodel fields
import numpy as np

##file i/o
def read_field(filename, nkx, nky, nz):
  data = np.fromfile(filename, dtype=float)
  data_real = np.reshape(data[0:nkx*nky*nz], (nky*nz, nkx))
  data_imag = np.reshape(data[nkx*nky*nz:2*nkx*nky*nz], (nky*nz, nkx))
  field = np.zeros((nkx, nky, nz), dtype=complex)
  for z in range(nz):
    field[:, :, z].real = data_real[z*nky:(z+1)*nky, :].T
    field[:, :, z].imag = data_imag[z*nky:(z+1)*nky, :].T
  return field

def write_field(filename, field):
  nkx, nky, nz = field.shape
  data_real = np.zeros((nky*nz, nkx))
  data_imag = np.zeros((nky*nz, nkx))
  for z in range(nz):
    data_real[z*nky:(z+1)*nky, :] = field[:, :, z].real.T
    data_imag[z*nky:(z+1)*nky, :] = field[:, :, z].imag.T
  data = np.zeros((2*nkx*nky*nz,))
  data[0:nkx*nky*nz] = np.reshape(data_real, (nkx*nky*nz,))
  data[nkx*nky*nz:2*nkx*nky*nz] = np.reshape(data_imag, (nkx*nky*nz,))
  f = open(filename, 'wb')  ##overwrites old file if exists
  data.tofile(f)
  f.close()

def output_ens(filename, xens):
  from netCDF4 import Dataset
  import os
  nens, ni, nj, nk = xens.shape
  if os.path.exists(filename):
    os.remove(filename)
  f = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
  ii = f.createDimension('i', ni)
  jj = f.createDimension('j', nj)
  kk = f.createDimension('k', nk)
  mm = f.createDimension('m', nens)
  dat = f.createVariable('xens', np.float32, ('m', 'k', 'j', 'i'))
  for n in range(nens):
    dat[n, :, :, :] = xens[n, :, :, :].T
  f.close()

##conversion between physical and spectral spaces
def fullspec(hfield):
  nkx, nky = hfield.shape
  kmax = nky-1
  hres = nkx+1
  sfield = np.zeros((hres, hres), dtype=complex)
  fup = hfield
  fup[kmax-1::-1, 0] = np.conjugate(fup[kmax+1:nkx, 0])
  fdn = np.conjugate(fup[::-1, nky-1:0:-1])
  sfield[1:hres, nky:hres] = fup
  sfield[1:hres, 1:nky] = fdn
  return sfield

def halfspec(sfield):
  n1, n2 = sfield.shape
  kmax = int(n1/2)-1
  hfield = sfield[1:n1, kmax+1:n2]
  return hfield

def spec2grid(fieldk):
  nkx, nky, nz = fieldk.shape
  nx = nkx+1
  ny = 2*nky
  field = np.zeros((nx, ny, nz))
  for z in range(nz):
    tmp = np.fft.ifftshift(fullspec(fieldk[:, :, z]))
    field[:, :, z] = nx*ny*np.real(np.fft.ifft2(tmp))
  return field

def grid2spec(field):
  nx, ny, nz = field.shape
  nkx = nx-1
  nky = int(ny/2)
  fieldk = np.zeros((nkx, nky, nz), dtype=complex)
  for z in range(nz):
    tmp = np.fft.fft2(field[:, :, z])/nx/ny
    fieldk[:, :, z] = halfspec(np.fft.fftshift(tmp))
  return fieldk

###spectral operation
def pwrspec2d(field):
  nx, ny, nz = field.shape
  FT = np.zeros((nx, ny, nz), dtype=complex)
  nupx = int(np.ceil((nx+1)/2))
  nupy = int(np.ceil((ny+1)/2))
  nup = max(nupx, nupy)
  wnx = generate_fft_index(nx)
  wny = generate_fft_index(ny)
  ky, kx = np.meshgrid(wnx, wny)
  k2d = np.sqrt((kx*(nup/nupx))**2 + (ky*(nup/nupy))**2)
  for z in range(nz):
    FT[:, :, z] = np.fft.fft2(field[:, :, z])
  P = (np.abs(FT)/nx/ny)**2
  wn = np.arange(0.0, nup)
  pwr = np.zeros((nup, nz))
  for z in range(nz):
    Pz = P[:, :, z]
    for w in range(nup):
      pwr[w, z] = np.sum(Pz[np.where(np.ceil(k2d)==w)])
  return wn, pwr

def spec_bandpass(xk, krange, s):
  kx_, ky_, z_ = get_coords(xk)
  Kh = np.sqrt(kx_**2 + ky_**2)
  xkout = xk.copy()
  r = scale_response(Kh, krange, s)
  return xkout * r

def scale_response(Kh, krange, s):
  ns = len(krange)
  r = np.zeros(Kh.shape)
  center_k = krange[s]
  if s == 0:
    r[np.where(Kh<=center_k)] = 1.0
  else:
    left_k = krange[s-1]
    ind = np.where(np.logical_and(Kh>=left_k, Kh<=center_k))
    r[ind] = np.cos((Kh[ind] - center_k)*(0.5*np.pi/(left_k - center_k)))**2
  if s == ns-1:
    r[np.where(Kh>=center_k)] = 1.0
  else:
    right_k = krange[s+1]
    ind = np.where(np.logical_and(Kh>=center_k, Kh<=right_k))
    r[ind] = np.cos((Kh[ind] - center_k)*(0.5*np.pi/(right_k - center_k)))**2
  return r

def spec_bandpass1(xk, Kmin, Kmax):
  kx_, ky_, z_ = get_coords(xk)
  Kh = np.sqrt(kx_**2 + ky_**2)
  xkout = xk.copy()
  xkout[np.where(Kh>Kmax)] = 0.0
  xkout[np.where(Kh<=Kmin)] = 0.0
  return xkout

##conversion between variables
def get_coords(psik):
  nkx, nky, nz = psik.shape
  kmax = nky-1
  kx_, ky_, z_ = np.mgrid[-kmax:kmax+1, 0:kmax+1, 0:nz]
  return kx_, ky_, z_

def psi2zeta(psik):
  kx_, ky_, z_ = get_coords(psik)
  zetak = -(kx_**2 + ky_**2) * psik
  return zetak

def psi2psi(psik):
  return psik

def psi2temp(psik):
  kx_, ky_, z_ = get_coords(psik)
  tempk = -np.sqrt(kx_**2 + ky_**2) * psik
  return tempk

def psi2u(psik):
  kx_, ky_, z_ = get_coords(psik)
  uk = -1j * ky_ * psik
  return uk

def psi2v(psik):
  kx_, ky_, z_ = get_coords(psik)
  vk = 1j * kx_ * psik
  return vk

def uv2zeta(uk, vk):
  kx_, ky_, z_ = get_coords(uk)
  zetak = 1j*kx_*vk - 1j*ky_*uk
  return zetak

def zeta2psi(zetak):
  nkx, nky, nz = zetak.shape
  kmax = nky-1
  kx_, ky_, z_ = get_coords(zetak)
  k2_ = kx_**2 + ky_**2
  k2_[kmax, 0, :] = 1  #set irrelavent point to 1 to avoid singularity in inversion
  psik = -(1.0/k2_) * zetak
  return psik

def temp2psi(tempk):
  nkx, nky, nz = tempk.shape
  kmax = nky-1
  kx_, ky_, z_ = get_coords(tempk)
  k1_ = np.sqrt(kx_**2 + ky_**2)
  k1_[kmax, 0, :] = 1
  psik = -(1.0/k1_) * tempk
  return psik

#another way to compute vorticity from psi
def vort(x):
  ni, nj = x.shape
  di = 2*np.pi/ni
  x11 = np.roll(np.roll(x, -1, axis=0), -1, axis=1)
  x12 = np.roll(np.roll(x, 1, axis=0), -1, axis=1)
  x21 = np.roll(np.roll(x, -1, axis=0), 1, axis=1)
  x22 = np.roll(np.roll(x, 1, axis=0), 1, axis=1)
  y = (x11 + x12 + x21 + x22 - 4*x) / di**2 /2
  return y

###spatial operation
def deriv_x(f):
  fx = 0.5*(np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0))
  return fx

def deriv_y(f):
  fy = 0.5*(np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1))
  return fy

def deriv_xy(f):
  fxy = 0.25*(np.roll(np.roll(f, -1, axis=0), -1, axis=1) + np.roll(np.roll(f, 1, axis=0), 1, axis=1) - np.roll(np.roll(f, -1, axis=0), 1, axis=1) - np.roll(np.roll(f, 1, axis=0), -1, axis=1))
  return fxy

def laplacian(f):
  del2f = (np.roll(f, -1, axis=0) + np.roll(f, 1, axis=0) + np.roll(f, -1, axis=1) + np.roll(f, 1, axis=1))/6 + (np.roll(np.roll(f, -1, axis=1), -1, axis=0) + np.roll(np.roll(f, -1, axis=1), 1, axis=0) + np.roll(np.roll(f, 1, axis=1), -1, axis=0) + np.roll(np.roll(f, 1, axis=1), 1, axis=0))/12 - f
  return del2f

def deriv_xx(f):
  fxx = (np.roll(f, -1, axis=0) + np.roll(f, 1, axis=0))/2 - f
  return fxx

def deriv_yy(f):
  fyy = (np.roll(f, -1, axis=1) + np.roll(f, 1, axis=1))/2 - f
  return fyy


def warp(Im, u, v):
  warp_Im = Im.copy()
  ni, nj = Im.shape
  for i in range(ni):
    for j in range(nj):
      warp_Im[i, j] = interp2d(Im, (i+u[i, j], j+v[i, j]))
  return warp_Im

def coarsen(Im, level):
  for k in range(level):
    ni, nj = Im.shape
    Im1 = 0.25*(Im[0:ni:2, :][:, 0:nj:2] + Im[1:ni:2, :][:, 0:nj:2] + Im[0:ni:2, 1:nj:2] + Im[1:ni:2, 1:nj:2])
    Im = Im1
  return Im

def sharpen(Im, level):
  for k in range(level):
    ni, nj = Im.shape
    Im1 = np.zeros((ni*2, nj))
    Im1[0:ni*2:2, :] = Im
    Im1[1:ni*2:2, :] = 0.5*(np.roll(Im, -1, axis=0) + Im)
    Im2 = np.zeros((ni*2, nj*2))
    Im2[:, 0:nj*2:2] = Im1
    Im2[:, 1:nj*2:2] = 0.5*(np.roll(Im1, -1, axis=1) + Im1)
    Im = Im2
  return Im

def interp2d(x, loc):
  ni, nj = x.shape
  io = loc[0]
  jo = loc[1]
  io1 = int(np.floor(io)) % ni
  jo1 = int(np.floor(jo)) % nj
  io2 = int(np.floor(io+1)) % ni
  jo2 = int(np.floor(jo+1)) % nj
  di = io - np.floor(io)
  dj = jo - np.floor(jo)
  xo = (1-di)*(1-dj)*x[io1, jo1] + di*(1-dj)*x[io2, jo1] + (1-di)*dj*x[io1, jo2] + di*dj*x[io2, jo2]
  return xo

def interp3d(x, loc):
  ni, nj, nk = x.shape
  io = loc[0]
  jo = loc[1]
  ko = loc[2]
  io1 = int(np.floor(io)) % ni
  jo1 = int(np.floor(jo)) % nj
  ko1 = int(np.floor(ko))
  io2 = int(np.floor(io+1)) % ni
  jo2 = int(np.floor(jo+1)) % nj
  ko2 = int(np.floor(ko+1))
  di = io - np.floor(io)
  dj = jo - np.floor(jo)
  dk = ko - np.floor(ko)
  xo1 = (1-di)*(1-dj)*x[io1, jo1, ko1] + di*(1-dj)*x[io2, jo1, ko1] + (1-di)*dj*x[io1, jo2, ko1] + di*dj*x[io2, jo2, ko1]
  xo2 = (1-di)*(1-dj)*x[io1, jo1, ko2] + di*(1-dj)*x[io2, jo1, ko2] + (1-di)*dj*x[io1, jo2, ko2] + di*dj*x[io2, jo2, ko2]
  xo = (1-dk)*xo1 + dk*xo2
  return xo

def regrid(x1, ni, nj):
  ni1, nj1 = x1.shape
  ii1, jj1 = np.mgrid[0:ni1:1.0*ni1/ni, 0:nj1:1.0*nj1/nj]
  x = np.zeros((ni, nj))
  for i in range(ni):
    for j in range(nj):
      x[i, j] = interp2d(x1, np.array([ii1[i, j], jj1[i, j]]))
  return x

def smooth(x, smth):
  if smth > 0:
    x_smooth = np.zeros(x.shape)
    cw = 0.0
    for i in np.arange(-smth, smth, 1):
      for j in np.arange(-smth, smth, 1):
        w = np.exp(-(i**2+j**2)/(smth/2.0)**2)
        cw += w
        x_smooth += w * np.roll(np.roll(x, j, axis=0), i, axis=1)
    x_smooth = x_smooth/cw
  else:
    x_smooth = x
  return x_smooth

def smooth_spec(wn, pwr, smth):
  n = wn.size
  pwr_smth = pwr.copy()
  for m in range(1, n):
    pwr_smth[m] = np.mean(pwr[int(max(0, np.floor(m/smth))):int(min(np.ceil(m*smth), n-1))+1])
  return pwr_smth

####random fields
def generate_fft_index(n):
  nup = int(np.ceil((n+1)/2))
  if n%2 == 0:
    wn = np.concatenate((np.arange(0, nup), np.arange(2-nup, 0)))
  else:
    wn = np.concatenate((np.arange(0, nup), np.arange(1-nup, 0)))
  return wn

def gaussian_random_field(Pk, n):
  wn = generate_fft_index(n)
  kx, ky = np.meshgrid(wn, wn)
  k2d = np.sqrt(kx**2 + ky**2)
  k2d[np.where(k2d==0.0)] = 1e-10
  noise = np.fft.fft2(np.random.normal(0, 1, (n, n)))
  amplitude = Pk(k2d)
  amplitude[np.where(k2d==1e-10)] = 0.0
  noise1 = np.real(np.fft.ifft2(noise * amplitude))
  return (noise1 - np.mean(noise1))/np.std(noise1)

####diagnostics
def rmse(x, xt):
  return np.sqrt(np.mean((x-xt)**2))

def sprd(xens):
  return np.sqrt(np.mean(np.std(xens, axis=0)**2))

def pattern_correlation(x, xt):
  nx, ny = x.shape
  x_mean = np.mean(x)
  xt_mean = np.mean(xt)
  xp = x - x_mean
  xtp = xt - xt_mean
  cov = np.sum(xp * xtp)
  x_norm = np.sum(xp ** 2)
  xt_norm = np.sum(xtp ** 2)
  pcorr = cov/np.sqrt(x_norm * xt_norm)
  return pcorr
