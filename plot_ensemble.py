#!/usr/bin/env python3
###read in ensemble state from bin files and plot spaghetti contours
import numpy as np
import util
import matplotlib.pyplot as plt
import matplotlib
import sys

if len(sys.argv) != 8:
  print("usage: ./plot_ensemble.py workdir casename kmax nz nens t ens_type")
  exit()

workdir = sys.argv[1]         #working directory
casename = sys.argv[2]        #casename for ensemble filtering experiment
kmax = int(sys.argv[3])       #maximum wavenumber (size of model n = 2*(kmax+1))
nz = int(sys.argv[4])         #number of vertical levels
nens = int(sys.argv[5])       #ensemble size
t = int(sys.argv[6])          #time step
ens_type = int(sys.argv[7])   #1: prior ensemble, 2: posterior ensemble

nkx = 2*kmax+1
nky = kmax+1
nx = 2*(kmax+1)
ny = 2*(kmax+1)
ii, jj = np.mgrid[0:nx, 0:ny]

lv = 0  #vertical level to plot
smth = 2  #smoothing the field for plotting
krange = (3, 8, 20)
s = 2

varname = r'$\theta$'
convertor = util.psi2temp
clevel_highlight = (-10, 10)
varmin = -40
varmax = 40
varint = 2

plt.switch_backend('Agg')
fig, ax = plt.subplots(2, 1, figsize=(5, 8))
matplotlib.rcParams['mathtext.fontset'] = 'cm'

def set_axis(ax, title):
  ax.set_aspect('equal', 'box')
  ax.set_xlim(0, nx)
  ax.set_ylim(0, ny)
  ax.set_xticks(np.arange(0, nx+1, 20))
  ax.set_yticks(np.arange(0, ny+1, 20))
  # ax.set_xlim(20, 60)
  # ax.set_ylim(40, 80)
  # ax.set_xticks(np.arange(20, 61, 10))
  # ax.set_yticks(np.arange(40, 81, 10))
  ax.tick_params(labelsize=10)
  ax.set_title(title, fontsize=20, fontname='cmr10')

def set_clevel(vout, vmin, vmax, vint):
  clevel = np.arange(vmin, vmax+vint, vint)
  vout[np.where(vout > vmax)] = vmax
  vout[np.where(vout < vmin)] = vmin
  return vout, clevel

def colorbar_ax(ax):
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  fig = ax.figure
  div = make_axes_locatable(ax)
  cax = div.append_axes('right', size='5%', pad=0.05)
  cax.tick_params(labelsize=10)
  return cax

###truth
filename = workdir+'/truth/{:05d}.bin'.format(t)
psik = util.read_field(filename, nkx, nky, nz)
# var = util.spec2grid(util.spec_bandpass(convertor(psik), krange, s))
var = util.spec2grid(convertor(psik))
var = np.roll(np.roll(var, -40, axis=0), 60, axis=1) #shift domain position for better plotting

out1 = util.smooth(var[:, :, lv], smth)
out1, clevel = set_clevel(out1, varmin, varmax, varint)
c = ax[0].contourf(ii, jj, out1, clevel, cmap='seismic')
ax[0].contour(ii, jj, out1, clevel_highlight, colors='black', linestyles='solid', linewidths=2)
cax = colorbar_ax(ax[0])
plt.colorbar(c, cax=cax)
set_axis(ax[0], varname+' truth')

###ensemble
if ens_type == 1:
  name = 'f_{:05d}'.format(t)
if ens_type == 2:
  name = '{:05d}'.format(t)
if ens_type == 3:
  name = 'fa_{:05d}'.format(t)

varens = np.zeros((nens, nx, ny, nz))
for m in range(nens):
  filename = workdir+'/'+casename+'/'+'{:04d}'.format(m+1)+'/'+name+'.bin'
  psik = util.read_field(filename, nkx, nky, nz)
  # varens[m, :, :, :] = util.spec2grid(util.spec_bandpass(convertor(psik), krange, s))
  varens[m, :, :, :] = util.spec2grid(convertor(psik))
varens = np.roll(np.roll(varens, -40, axis=1), 60, axis=2) #shift domain position for better plotting

##some error statistics
# for m in range(nens):
#   varmean = np.mean(varens[m:m+1, :, :, lv], axis=0)
#   rmse = util.rmse(varmean, var[:, :, lv])
#   pcorr = util.pattern_correlation(varmean, var[:, :, lv])
#   print('error = {:7.2f}'.format(rmse))
#   print('pattern correlation = {:7.2f}'.format(pcorr))
varmean = np.mean(varens[:, :, :, lv], axis=0)
rmse = util.rmse(varmean, var[:, :, lv])
pcorr = util.pattern_correlation(varmean, var[:, :, lv])
print('error = {:7.2f}'.format(rmse))
print('pattern correlation = {:7.2f}'.format(pcorr))
sprd = util.sprd(varens[:, :, :, lv])
print('ensemble spread = {:7.2f}'.format(sprd))


cmap = [plt.cm.jet(m) for m in np.linspace(0, 1, nens)]
for m in range(nens):
  out = util.smooth(varens[m, :, :, lv], smth)
  ax[1].contour(ii, jj, out, clevel_highlight, colors=[cmap[m][0:3]], linestyles='solid', linewidths=1)
ax[1].contour(ii, jj, out1, clevel_highlight, colors='black', linestyles='solid', linewidths=2)
cax = colorbar_ax(ax[1])
cax.set_visible(False)
set_axis(ax[1], varname+' ensemble')

fig.tight_layout()
# plt.savefig(workdir+'/'+casename+'/ensemble_'+name+'.pdf')
plt.savefig('1.pdf')
