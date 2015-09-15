# update PDF of ds and icorr
import sys
import os
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
plt.ion()

from constants import *
from functions.function import weightedAvgAndStd
from pyre.distributions import Lognormal

#NO_EVIDENCE_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'no_evidence')
#EVIDENCE2_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_condition_state' )

NO_EVIDENCE_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'no_evidence_shear')
EVIDENCE2_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_condition_state_shear' )
# NO_EVIDENCE_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'no_evidence')
# EVIDENCE2_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_condition_state' )

EVIDENCE3_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_half_cell_shear')
EVIDENCE4_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_corrosion_rate_shear' )

#if __name__ == '__main__':
datafile = os.path.join(NO_EVIDENCE_DATA_PATH, 'LWS_results.txt')
service_time = np.loadtxt(datafile)[0,:]
    
display_yr = 78
    
# corrosion rate data (since icorr has model error, this model error should be divided that shouldn't be updated)
icorr_no_evidence = np.load( os.path.join(NO_EVIDENCE_DATA_PATH, 'corrosion_rate_history.npy') )[service_time==display_yr, :].flatten()
icorr_evidence2 = np.load( os.path.join(EVIDENCE2_DATA_PATH, 'corrosion_rate_history.npy') )[service_time==display_yr, :].flatten()
icorr_evidence3 = np.load( os.path.join(EVIDENCE3_DATA_PATH, 'corrosion_rate_history.npy') )[service_time==display_yr, :].flatten()
icorr_evidence4 = np.load( os.path.join(EVIDENCE4_DATA_PATH, 'corrosion_rate_history.npy') )[service_time==display_yr, :].flatten()
    
# residual diameter data
ds_no_evidence = np.load( os.path.join(NO_EVIDENCE_DATA_PATH, 'residual_diameter_history.npy') )[service_time==display_yr, :].flatten()
ds_evidence2 = np.load( os.path.join(EVIDENCE2_DATA_PATH, 'residual_diameter_history.npy') )[service_time==display_yr, :].flatten()
ds_evidence3 = np.load( os.path.join(EVIDENCE3_DATA_PATH, 'residual_diameter_history.npy') )[service_time==display_yr, :].flatten()
ds_evidence4 = np.load( os.path.join(EVIDENCE4_DATA_PATH, 'residual_diameter_history.npy') )[service_time==display_yr, :].flatten()

# residual flexural strength data
rc_flex_no_evidence = np.load( os.path.join(NO_EVIDENCE_DATA_PATH, 'rc_flexure_history.npy') )[service_time==display_yr, :].flatten()
rc_flex_evidence2 = np.load( os.path.join(EVIDENCE2_DATA_PATH, 'rc_flexure_history.npy') )[service_time==display_yr, :].flatten()
rc_flex_evidence3 = np.load( os.path.join(EVIDENCE3_DATA_PATH, 'rc_flexure_history.npy') )[service_time==display_yr, :].flatten()
rc_flex_evidence4 = np.load( os.path.join(EVIDENCE4_DATA_PATH, 'rc_flexure_history.npy') )[service_time==display_yr, :].flatten()

# residual shear strength data
rc_shear_no_evidence = np.load( os.path.join(NO_EVIDENCE_DATA_PATH, 'rc_shear_history.npy') )[service_time==display_yr, :].flatten()
rc_shear_evidence2 = np.load( os.path.join(EVIDENCE2_DATA_PATH, 'rc_shear_history.npy') )[service_time==display_yr, :].flatten()
rc_shear_evidence3 = np.load( os.path.join(EVIDENCE3_DATA_PATH, 'rc_shear_history.npy') )[service_time==display_yr, :].flatten()
rc_shear_evidence4 = np.load( os.path.join(EVIDENCE4_DATA_PATH, 'rc_shear_history.npy') )[service_time==display_yr, :].flatten()
    
# weight data
weights1 = np.ones(np.shape(ds_no_evidence))
weights2 = np.load( os.path.join(EVIDENCE2_DATA_PATH, 'likelihood_weighting.npy') )
weights3 = np.load( os.path.join(EVIDENCE3_DATA_PATH, 'likelihood_weighting.npy') )
weights4 = np.load( os.path.join(EVIDENCE4_DATA_PATH, 'likelihood_weighting.npy') )

# fitting data with lognormal distribution
# for ds
ds_mu1, ds_sigma1 = weightedAvgAndStd(ds_no_evidence, weights1)
ds_mu2, ds_sigma2 = weightedAvgAndStd(ds_evidence2, weights2)
ds_mu3, ds_sigma3 = weightedAvgAndStd(ds_evidence3, weights3)
ds_mu4, ds_sigma4 = weightedAvgAndStd(ds_evidence4, weights4)
ds_lognormal1 = Lognormal('no_ev', ds_mu1, ds_sigma1)
ds_lognormal2 = Lognormal('ev2', ds_mu2, ds_sigma2)
ds_lognormal3 = Lognormal('ev3', ds_mu3, ds_sigma3)
ds_lognormal4 = Lognormal('ev4', ds_mu4, ds_sigma4)
# for flexure resistance
flex_mu1, flex_sigma1 = weightedAvgAndStd(rc_flex_no_evidence, weights1)
flex_mu2, flex_sigma2 = weightedAvgAndStd(rc_flex_evidence2, weights2)
flex_mu3, flex_sigma3 = weightedAvgAndStd(rc_flex_evidence3, weights3)
flex_mu4, flex_sigma4 = weightedAvgAndStd(rc_flex_evidence4, weights4)
flex_lognormal1 = Lognormal('no_ev', flex_mu1, flex_sigma1)
flex_lognormal2 = Lognormal('ev2', flex_mu2, flex_sigma2)
flex_lognormal3 = Lognormal('ev3', flex_mu3, flex_sigma3)
flex_lognormal4 = Lognormal('ev4', flex_mu4, flex_sigma4)
# for shear resistance
shear_mu1, shear_sigma1 = weightedAvgAndStd(rc_shear_no_evidence, weights1)
shear_mu2, shear_sigma2 = weightedAvgAndStd(rc_shear_evidence2, weights2)
shear_mu3, shear_sigma3 = weightedAvgAndStd(rc_shear_evidence3, weights3)
shear_mu4, shear_sigma4 = weightedAvgAndStd(rc_shear_evidence4, weights4)
shear_lognormal1 = Lognormal('no_ev', shear_mu1, shear_sigma1)
shear_lognormal2 = Lognormal('ev2', shear_mu2, shear_sigma2)
shear_lognormal3 = Lognormal('ev3', shear_mu3, shear_sigma3)
shear_lognormal4 = Lognormal('ev4', shear_mu4, shear_sigma4)

## =======================================================================##
## generate and save figures
## =======================================================================##
    
plt.close('all')
plt.rc('font', family='serif', size=12)
all_figures = []
num_bins = 50
    
# results of icorr
f = plt.figure()
plt.hold(True)
# histogram: icorr
dummy, bins_evidence, dummy = plt.hist(icorr_no_evidence[icorr_no_evidence!=0], bins=num_bins, normed=1, color='blue', histtype='step')
dummy, bins_evidence2, dummy = plt.hist(icorr_evidence2[icorr_evidence2!=0], bins=num_bins, weights=weights2[icorr_evidence2!=0], normed=1, color='red', histtype='step')
dummy, bins_evidence3, dummy = plt.hist(icorr_evidence3[icorr_evidence3!=0], bins=num_bins, weights=weights3[icorr_evidence3!=0], normed=1, color='green', histtype='step')
dummy, bins_evidence4, dummy = plt.hist(icorr_evidence4[icorr_evidence4!=0], bins=num_bins, weights=weights4[icorr_evidence4!=0], normed=1, color='k', histtype='step')
plt.xlabel(r'corrosion rate ($\mu A/cm^2$)')
plt.ylabel(r'PDF')
plt.grid(True)
plt.legend(loc='upper left', prop={'size':12})
all_figures.append(f)

# histogram: residual diameter
f = plt.figure()
plt.hold(True)
# histogram: ds
dummy, bins_evidence1, dummy = plt.hist(ds_no_evidence[ds_no_evidence!=DS_REGION1_MEAN], bins=num_bins*3, normed=1, color='blue', histtype='step')
dummy, bins_evidence2, dummy = plt.hist(ds_evidence2[ds_evidence2!=DS_REGION1_MEAN], bins=num_bins, weights=weights2[ds_evidence2!=DS_REGION1_MEAN], normed=1, color='red', histtype='step')
dummy, bins_evidence3, dummy = plt.hist(ds_evidence3[ds_evidence3!=DS_REGION1_MEAN], bins=num_bins*3, weights=weights3[ds_evidence3!=DS_REGION1_MEAN], normed=1, color='green', histtype='step')
dummy, bins_evidence4, dummy = plt.hist(ds_evidence4[ds_evidence4!=DS_REGION1_MEAN], bins=num_bins*3, weights=weights4[ds_evidence4!=DS_REGION1_MEAN], normed=1, color='yellow', histtype='step')
#ds_str1 = np.linspace(bins_evidence1[0], bins_evidence1[-1], 1000)
#ds_str2 = np.linspace(bins_evidence2[0], bins_evidence2[-1], 1000)
#ds_str3 = np.linspace(bins_evidence3[0], bins_evidence3[-1], 1000)
#ds_str4 = np.linspace(bins_evidence4[0], bins_evidence4[-1], 1000)
#plt.plot(ds_str1, ds_lognormal1.rv.pdf(ds_str1), 'b-', lw=1.5)
#plt.plot(ds_str2, ds_lognormal2.rv.pdf(ds_str2), 'r-', lw=1.5)
#plt.plot(ds_str3, ds_lognormal3.rv.pdf(ds_str3), 'g-', lw=1.5)
#plt.plot(ds_str4, ds_lognormal4.rv.pdf(ds_str4), 'y-', lw=1.5)
plt.xlim( (6, 13) )
#plt.ylim( (0, 3) )
ax = plt.gca()
ax.annotate('Without evidence (year 78)', xy=(11.5, 0.6), xytext=(8.2, 0.9), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
ax.annotate('Condition states (year 78)', xy=(10., 0.4), xytext=(6.7, 0.7), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
ax.annotate('Half-cell (year 78)', xy=(12.0, 1.2), xytext=(9.4, 1.5), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
ax.annotate('Corrosion rate (year 78)', xy=(11.3, 1.0), xytext=(8.2, 1.3), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))

plt.xlabel(r'residual diameter of corroded rebars ($mm$)')
plt.ylabel(r'PDF')
plt.grid(False)
plt.legend(loc='upper left', prop={'size':12})
all_figures.append(f)


# histogram: residual flexural strength
f = plt.figure()
plt.hold(True)
# histogram: mu
dummy, bins_evidence1, dummy = plt.hist(rc_flex_no_evidence, bins=num_bins, normed=1, edgecolor='blue', color='white')
dummy, bins_evidence2, dummy = plt.hist(rc_flex_evidence2, bins=num_bins, weights=weights2, normed=1, edgecolor='red', color='white')
#dummy, bins_evidence3, dummy = plt.hist(rc_flex_evidence3, bins=num_bins, weights=weights3, normed=1, color='green', histtype='step')
#dummy, bins_evidence4, dummy = plt.hist(rc_flex_evidence4, bins=num_bins, weights=weights4, normed=1, color='k', histtype='step')
flex_str1 = np.linspace(bins_evidence1[0], bins_evidence1[-1], 1000)
flex_str2 = np.linspace(bins_evidence2[0], bins_evidence2[-1], 1000)
#flex_str3 = np.linspace(bins_evidence3[0], bins_evidence3[-1], 1000)
#flex_str4 = np.linspace(bins_evidence4[0], bins_evidence4[-1], 1000)
plt.plot(flex_str1, flex_lognormal1.rv.pdf(flex_str1), 'b-', lw=1.5)
plt.plot(flex_str2, flex_lognormal2.rv.pdf(flex_str2), 'r-', lw=1.5)
#plt.plot(flex_str3, flex_lognormal3.rv.pdf(flex_str3), 'g-', lw=1.5)
#plt.plot(flex_str4, flex_lognormal4.rv.pdf(flex_str4), 'y-', lw=1.5)
ax = plt.gca()
ax.annotate('condition states\n(year 78)', xy=(1331.00, 1.45e-3), xytext=(600, 1.71e-3),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
ax.annotate('w/o evidence (year 78)', xy=(1693.5, 1.51e-3), xytext=(1953.6, 1.72e-3),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
LogFit = plt.Line2D((0,1),(0,0), color='k', ls='-')
plt.legend([LogFit], [r'Lognormal fitting'], fontsize='12', bbox_to_anchor=(1, 1), frameon=False)
plt.xlabel(r'flexural capacity ($kN \cdot m$)')
plt.ylabel(r'PDF')
plt.grid(False)
plt.legend(loc='upper left', prop={'size':12})
all_figures.append(f)


# histogram: residual shear strength
f = plt.figure()
plt.hold(True)
# histogram: icorr
dummy, bins_evidence1, dummy = plt.hist(rc_shear_no_evidence, bins=num_bins, normed=1, edgecolor='blue', color='white')
dummy, bins_evidence2, dummy = plt.hist(rc_shear_evidence2, bins=num_bins, weights=weights2, normed=1, edgecolor='red', color='white')
#dummy, bins_evidence3, dummy = plt.hist(rc_shear_evidence3, bins=num_bins, weights=weights3, normed=1, edgecolor='green', color='white')
#dummy, bins_evidence4, dummy = plt.hist(rc_shear_evidence4, bins=num_bins, weights=weights4, normed=1, edgecolor='yellow', color='white')
shear_str1 = np.linspace(bins_evidence1[0], bins_evidence1[-1], 1000)
shear_str2 = np.linspace(bins_evidence2[0], bins_evidence2[-1], 1000)
#shear_str3 = np.linspace(bins_evidence3[0], bins_evidence3[-1], 1000)
#shear_str4 = np.linspace(bins_evidence4[0], bins_evidence4[-1], 1000)
plt.plot(shear_str1, shear_lognormal1.rv.pdf(shear_str1), 'b-', lw=1.5)
plt.plot(shear_str2, shear_lognormal2.rv.pdf(shear_str2), 'r-', lw=1.5)
#plt.plot(shear_str3, shear_lognormal3.rv.pdf(shear_str3), 'g-', lw=1.5)
#plt.plot(shear_str4, shear_lognormal4.rv.pdf(shear_str4), 'y-', lw=1.5)

ax = plt.gca()
ax.annotate('w/o evidence (year 78)', xy=(793, 0.00289), xytext=(1043.48, 3.43e-3),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
ax.annotate('condition states (year 78)', xy=(660.81, 0.0036), xytext=(911.29, 0.00414),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
LogFit = plt.Line2D((0,1),(0,0), color='k', ls='-')
plt.legend([LogFit], [r'Lognormal fitting'], fontsize='12', bbox_to_anchor=(1, 0.1), frameon=False)

plt.xlabel(r'shear capacity ($kN$)')
plt.ylabel(r'PDF')
plt.grid(False)
plt.legend(loc='upper left', prop={'size':12})
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
