# update PDF of ds and icorr
import sys
import os
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

from constants import *
from functions.function import weightedAvgAndStd
from pyre.distributions import Lognormal

NO_STR_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_condition_state_shear')
STR1_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'frp_U_anchor_yr78', 'evidence_condition_state_shear')
STR2_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'frp_U_bond_yr78', 'evidence_condition_state_shear')

#if __name__ == '__main__':
datafile = os.path.join(NO_STR_DATA_PATH, 'LWS_results.txt')
service_time = np.loadtxt(datafile)[0,:]
    
display_yr = 78 

# residual shear strength data
rc_shear_no_str = np.load( os.path.join(NO_STR_DATA_PATH, 'rc_shear_history.npy') )[service_time==display_yr, :].flatten()
rc_shear_str1 = np.load( os.path.join(STR1_DATA_PATH, 'rc_shear_history.npy') )[service_time==display_yr, :].flatten()
rc_shear_str2 = np.load( os.path.join(STR2_DATA_PATH, 'rc_shear_history.npy') )[service_time==display_yr, :].flatten()
    
# weight data
weights0 = np.load( os.path.join(NO_STR_DATA_PATH, 'likelihood_weighting.npy') )
weights1 = np.load( os.path.join(STR1_DATA_PATH, 'likelihood_weighting.npy') )
weights2 = np.load( os.path.join(STR2_DATA_PATH, 'likelihood_weighting.npy') )

## parameters and fitting distribution
mu0, sigma0 = weightedAvgAndStd(rc_shear_no_str, weights0)
mu1, sigma1 = weightedAvgAndStd(rc_shear_str1, weights1)
mu2, sigma2 = weightedAvgAndStd(rc_shear_str2, weights2)
lognormal0 = Lognormal('str0', mu0, sigma0)
lognormal1 = Lognormal('str1', mu1, sigma1)
lognormal2 = Lognormal('str2', mu2, sigma2)
print 'Before strengthening: mean = {}, std = {}'.format(mu0, sigma0)
print 'U-jacket with anchor: mean = {}, std = {}'.format(mu1, sigma1)
print 'U-jacket without anhor: mean = {}, std = {}'.format(mu2, sigma2)

## =======================================================================##
## generate and save figures
## =======================================================================##
    
plt.close('all')
plt.rc('font', family='serif', size=12)
all_figures = []
num_bins = 50

# histogram: residual shear strength
f = plt.figure()
plt.hold(True)
# histogram: icorr
#histtype='step'
dummy, bins_no_str, dummy = plt.hist(rc_shear_no_str, bins=num_bins, weights=weights0, normed=1, color='white', edgecolor='b', alpha=1)
dummy, bins_str1, dummy = plt.hist(rc_shear_str1, bins=num_bins, weights=weights1, normed=1, color='white', edgecolor='r', alpha=1)
dummy, bins_str2, dummy = plt.hist(rc_shear_str2, bins=num_bins, weights=weights2, normed=1, color='white', edgecolor='g', alpha=1)
shear_str0 = np.linspace(bins_no_str[0], bins_no_str[-1], 1000)
shear_str1 = np.linspace(bins_str1[0], bins_str1[-1], 1000)
shear_str2 = np.linspace(bins_str2[0], bins_str2[-1], 1000)
plt.plot(shear_str0, lognormal0.rv.pdf(shear_str0), 'b-', lw=1.5)
plt.plot(shear_str1, lognormal1.rv.pdf(shear_str1), 'r-', lw=1.5)
plt.plot(shear_str2, lognormal2.rv.pdf(shear_str2), 'g-', lw=1.5)

ax = plt.gca()
ax.annotate('Before strengthening', xy=(713.71, 0.003), xytext=(1173.4, 0.00375), 
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
ax.annotate('U-jacketing w/ anchors', xy=(955.645, 0.00173), xytext=(1415.335, 2.48e-3), 
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
ax.annotate('U-jacketing w/o anchors', xy=(1451.2, 5.15e-4), xytext=(1862.92, 1.254e-3), 
            arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
LogFit = plt.Line2D((0,1),(0,0), color='k', ls='-')
plt.legend([LogFit], [r'Lognormal fitting'], fontsize='12', bbox_to_anchor=(1, 1), frameon=False)
plt.xlabel(r'shear capacity ($kN$)')
plt.ylabel(r'PDF')
plt.grid(False)
#plt.legend(loc='upper left', prop={'size':12})
all_figures.append(f)

    
## show figures
plt.show()
    
# save figures
isSave = raw_input('Save figures? (y/n)')
if isSave == 'y':
    plt.close('all')
    fig_no = 1
    for fig in all_figures:
        fig.savefig('fig'+str(fig_no)+'.eps')
        fig_no += 1
    print 'figures saveed. Postprocessing ends.'
elif isSave == 'n':
    print 'Postprocessing ends.'
else:
    print 'Illegal input: figures not saved'
