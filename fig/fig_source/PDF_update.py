# update PDF of ds and icorr
import sys
import os
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

NO_EVIDENCE_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'no_evidence_shear')
EVIDENCE2_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_condition_state_shear' )
EVIDENCE3_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_half_cell_shear')
EVIDENCE4_DATA_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_corrosion_rate_shear' )

#if __name__ == '__main__':
datafile = os.path.join(NO_EVIDENCE_DATA_PATH, 'LWS_results.txt')    
service_time = np.loadtxt(datafile)[0,:]
    
display_yr = 80
    
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
weights2 = np.load( os.path.join(EVIDENCE2_DATA_PATH, 'likelihood_weighting.npy') )
weights3 = np.load( os.path.join(EVIDENCE3_DATA_PATH, 'likelihood_weighting.npy') )
weights4 = np.load( os.path.join(EVIDENCE4_DATA_PATH, 'likelihood_weighting.npy') )



## =======================================================================##
## generate and save figures
## =======================================================================##
    
plt.close('all')
plt.rc('font', family='serif', size=12)
all_figures = []
num_bins = 100
    
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

# converge of icorr (evidence 2)



# histogram: residual diameter
f = plt.figure()
plt.hold(True)
# histogram: icorr
dummy, bins_evidence, dummy = plt.hist(ds_no_evidence[ds_no_evidence!=16], bins=num_bins*5, normed=1, color='blue', histtype='step')
dummy, bins_evidence2, dummy = plt.hist(ds_evidence2[ds_evidence2!=16], bins=num_bins, weights=weights2[ds_evidence2!=16], normed=1, color='red', histtype='step')
dummy, bins_evidence3, dummy = plt.hist(ds_evidence3[ds_evidence3!=16], bins=num_bins*5, weights=weights3[ds_evidence3!=16], normed=1, color='green', histtype='step')
dummy, bins_evidence4, dummy = plt.hist(ds_evidence4[ds_evidence4!=16], bins=num_bins*5, weights=weights4[ds_evidence4!=16], normed=1, color='k', histtype='step')
plt.xlim( (13, 16) )
plt.ylim( (0, 3) )
plt.xlabel(r'residual diameter ($mm$)')
plt.ylabel(r'PDF')
plt.grid(True)
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