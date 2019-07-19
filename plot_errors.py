#!/usr/bin/env python3
import numpy as np
import util
import matplotlib.pyplot as plt
import matplotlib
import sys

workdir = '/glade/scratch/mying/qgmodel_enkf/multiscale'

###Case 1: changing ensemble size N=5-80
cp = np.array([2, 2, 2, 2, 2])
nens = np.array([5, 10, 20, 40, 80])
###Case 2: changing cycling period cp=1-10
# cp = np.array([1, 2, 3, 5, 8, 10])
# nens = np.array([20, 20, 20, 20, 20, 20])
nt = 200
dt = 0.05

def ttest(a, b):
  from scipy import stats
  N = a.size
  var_a = a.var(ddof=1)
  var_b = b.var(ddof=1)
  s = np.sqrt((var_a + var_b)/2)
  diff = a.mean() - b.mean()
  t = diff/(s * np.sqrt(2/N))
  df = 2*N - 2
  if diff < 0:
    p_value = stats.t.cdf(t, df=df)
  else:
    p_value = 1 - stats.t.cdf(t, df=df)
  return diff, p_value

rmse_sample = np.zeros((nt, cp.size, 2, 3))
sprd_sample = np.zeros((nt, cp.size, 2, 3))
for i in range(cp.size):
  rmse = np.load(workdir+'/cp{:d}_ss_N{:d}'.format(cp[i], nens[i])+'/rmse.npy')
  rmse_sample[:, i, :, 0] = rmse[0:nt, :]
  rmse = np.load(workdir+'/cp{:d}_ms_N{:d}'.format(cp[i], nens[i])+'/rmse.npy')
  rmse_sample[:, i, :, 1] = rmse[0:nt, :]
  rmse = np.load(workdir+'/cp{:d}_msa_N{:d}'.format(cp[i], nens[i])+'/rmse.npy')
  rmse_sample[:, i, :, 2] = rmse[0:nt, :]
  sprd = np.load(workdir+'/cp{:d}_ss_N{:d}'.format(cp[i], nens[i])+'/sprd.npy')
  sprd_sample[:, i, :, 0] = sprd[0:nt, :]
  sprd = np.load(workdir+'/cp{:d}_ms_N{:d}'.format(cp[i], nens[i])+'/sprd.npy')
  sprd_sample[:, i, :, 1] = sprd[0:nt, :]
  sprd = np.load(workdir+'/cp{:d}_msa_N{:d}'.format(cp[i], nens[i])+'/sprd.npy')
  sprd_sample[:, i, :, 2] = sprd[0:nt, :]

plt.switch_backend('Agg')
fig, ax = plt.subplots(1, 2, figsize=(12, 4))

for j in range(2):
  for i in range(cp.size):
    bx = ax[j].boxplot(rmse_sample[:, i, 1-j, 0], positions=[(i+1)*2-0.3], widths=0.2, showfliers=False, patch_artist=True)
    for item in ['boxes', 'whiskers', 'medians', 'caps']:
      plt.setp(bx[item], color='k')
    plt.setp(bx['boxes'], facecolor=(0.5, 0.5, 0.5))
    bx = ax[j].boxplot(rmse_sample[:, i, 1-j, 1], positions=[(i+1)*2+0.0], widths=0.2, showfliers=False, patch_artist=True)
    for item in ['boxes', 'whiskers', 'medians', 'caps']:
      plt.setp(bx[item], color='b')
    plt.setp(bx['boxes'], facecolor=(0.5, 0.5, 1.0))
    bx = ax[j].boxplot(rmse_sample[:, i, 1-j, 2], positions=[(i+1)*2+0.3], widths=0.2, showfliers=False, patch_artist=True)
    for item in ['boxes', 'whiskers', 'medians', 'caps']:
      plt.setp(bx[item], color='r')
    plt.setp(bx['boxes'], facecolor=(1.0, 0.5, 0.5))
  # ax[j].set_xlim(1, 13)
  # ax[j].set_xticks(np.arange(2, 13, 2))
  # ax[j].set_xticklabels(('0.05', '0.1', '0.15', '0.25', '0.4', '0.5'))
  ax[j].set_xlim(1, 11)
  ax[j].set_xticks(np.arange(2, 11, 2))
  ax[j].set_xticklabels(('5', '10', '20', '40', '80'))
  ax[j].tick_params(labelsize=14)
# ax[0].set_ylim(0.8, 2.2)
# ax[1].set_ylim(0.5, 5.5)
ax[0].set_ylim(0.8, 2.4)
ax[1].set_ylim(1.0, 3.0)

plt.savefig('1.pdf')

for j in range(2):
  if j == 0:
    print('Posterior')
  if j == 1:
    print('Prior')
  for i in range(cp.size):
    # print('  CP = {}'.format(cp[i]))
    print('  N = {}'.format(nens[i]))
    rmseout = np.mean(rmse_sample[:, i, 1-j, :], axis=0)
    print('    RMSE={}'.format(rmseout))
    sprdout = np.mean(sprd_sample[:, i, 1-j, :], axis=0)
    print('    CR  ={}'.format(sprdout/rmseout))
    diff, pval = ttest(rmse_sample[:, i, 1-j, 1], rmse_sample[:, i, 1-j, 0])
    print('    MS-SS = {:7.3f}, p = {:5.3f}'.format(diff, pval))
    diff, pval = ttest(rmse_sample[:, i, 1-j, 2], rmse_sample[:, i, 1-j, 1])
    print('    MSA-MS = {:7.3f}, p = {:5.3f}'.format(diff, pval))

