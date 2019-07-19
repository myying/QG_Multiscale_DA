#!/usr/bin/env python3
###read in ensemble state from bin files and plot error spectrum
import numpy as np
import util
import matplotlib.pyplot as plt
import sys

workdir = '/glade/scratch/mying/qgmodel_enkf/multiscale_dt0.05'
obsdir = 'obs'
ens_type = 'post'
# t1 = int(sys.argv[1])
# t2 = int(sys.argv[2])
t1 = 0
t2 = 999
tint = 5
cp = 10

power_law = -3.0
lv = 0  #vertical level to plot
smth = 1.3

plt.switch_backend('Agg')
ax = plt.subplot(111)

pwr = np.load(workdir+'/truth/spectrum.npy')

nup, nz, nt = pwr.shape
wn = np.arange(0.0, nup)
ax.semilogx(wn, util.smooth_spec(wn, np.mean(pwr[:, lv, t1:t2:tint], axis=1), smth), color='black', label='Reference')

obserr = np.load(workdir+'/'+obsdir+'/spectrum_obserr.npy')
nupob, nzob, ntob = obserr.shape
wnob = np.arange(0.0, nupob)
ax.semilogx(wnob, util.smooth_spec(wnob, np.mean(obserr[:, 0, t1:t2:tint], axis=1), smth), color=[.7, .7, .7], label='Obs Error')

casename = 'cp{:d}_disp'.format(cp)
pwr_ens = np.load(workdir+'/'+casename+'/spectrum_'+ens_type+'_ens.npy')
err = np.load(workdir+'/'+casename+'/spectrum_'+ens_type+'_err.npy')
ax.semilogx(wn, util.smooth_spec(wn, np.mean(err[:, lv, t1:t2:tint], axis=1), smth), color='r', label='MSA')

casename = 'cp{:d}_nodisp'.format(cp)
pwr_ens = np.load(workdir+'/'+casename+'/spectrum_'+ens_type+'_ens.npy')
err = np.load(workdir+'/'+casename+'/spectrum_'+ens_type+'_err.npy')
ax.semilogx(wn, util.smooth_spec(wn, np.mean(err[:, lv, t1:t2:tint], axis=1), smth), color='b', linestyle='-', label='MS')
# sprd = np.mean(pwr_ens, axis=2) #ensemble spread spectrum
# ax.semilogx(wn, util.smooth_spec(wn, np.mean(sprd[:, lv, t1:t2:tint], axis=1), smth), color='r', linestyle=':')
# casename = 'enkf_disp'
# pwr_ens = np.load(workdir+'/'+casename+'/spectrum_'+ens_type+'_ens.npy')
# err = np.load(workdir+'/'+casename+'/spectrum_'+ens_type+'_err.npy')
# ax.semilogx(wn, util.smooth_spec(wn, np.mean(err[:, lv, t1:t2:tint], axis=1), smth), color=[.3, .7, .3], linestyle='-', label='EnKF+Align')

# wn[0]=1e-10
# ax.semilogx(wn, 1e1*wn**power_law, color=[.7, .7, .7], label=r'$k^{'+'{:4.1f}'.format(power_law)+'}$') #reference line of power law

ax.grid(True, which='both', axis='x', linestyle=':', color=[.7, .7, .7])
ax.grid(True, which='major', axis='y', linestyle=':', color=[.7, .7, .7])
ax.set_xlim(1, 100)
ax.set_ylim(0, 0.2)
#ax.set_ylim(1e-4, 1e1)
ax.tick_params(labelsize=12)
ax.set_xlabel('wavenumber', fontsize=16)
ax.set_title(r'$\theta$ spectrum', fontsize=20)
ax.legend(loc=0, fontsize=12)

# plt.savefig(workdir+'/spectrum{:05d}.png'.format(t1))
plt.savefig('1.pdf')
