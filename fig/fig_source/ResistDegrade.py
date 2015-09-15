# check cl and cdf of corrosion initiation cdf
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

FIG_PATH = os.path.join(os.path.abspath('./'), 'fig')  
FIG_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data') 
    
if __name__ == '__main__':
    ## =======================================================================##
    ## get results without evidence
    ## =======================================================================##
    #datafile = FIG_DATAFILE_PATH+'no_evidence/LWS_results.txt'  
    datafile = os.path.join(FIG_DATAFILE_PATH, 'no_evidence_shear', 'LWS_results.txt')
    service_time = np.loadtxt(datafile)[0,:]
    
    yrOfInterest = 100.0
    end_index = np.where(service_time==yrOfInterest)[0]+1
    
    # residual diameter history
    rc_shear_mean_history = np.loadtxt(datafile)[21,:]
    rc_shear_std_history = np.loadtxt(datafile)[22,:]
    rc_shear_cov_history = rc_shear_std_history/rc_shear_mean_history
    datafile = os.path.join(FIG_DATAFILE_PATH, 'no_evidence', 'LWS_results.txt')
    rc_flex_mean_history = np.loadtxt(datafile)[19,:]
    rc_flex_std_history = np.loadtxt(datafile)[20,:]
    rc_flex_cov_history = rc_flex_std_history/rc_flex_mean_history

    
    ## =======================================================================##
    ## get results with evidence
    ## =======================================================================##

    # shear
    datafile = os.path.join(FIG_DATAFILE_PATH, 'evidence_condition_state_shear', 'LWS_results.txt')
    rc_shear_mean_condition_state= np.loadtxt(datafile)[21,:]
    rc_shear_std_condition_state = np.loadtxt(datafile)[22,:]
    rc_shear_cov_condition_state = rc_shear_std_condition_state/rc_shear_mean_condition_state
    # flexure
    datafile = os.path.join(FIG_DATAFILE_PATH, 'evidence_condition_state', 'LWS_results.txt')
    rc_flex_mean_condition_state = np.loadtxt(datafile)[19,:]
    rc_flex_std_condition_state = np.loadtxt(datafile)[20,:]
    rc_flex_cov_condition_state = rc_flex_std_condition_state/rc_flex_mean_condition_state

    ## =======================================================================##
    ## generate and save figures
    ## =======================================================================##
    
    plt.close('all')
    plt.rc('font', family='serif', size=12)
    all_figures = []
    
    # residual diameter mean
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    # MC reults
    ax1.plot(service_time[:end_index], rc_flex_mean_history[:end_index]/rc_flex_mean_history[0], 'b--', label='mean (w/o evidence)', linewidth=1.5 )
    ax1.plot(service_time[:end_index], rc_shear_mean_history[:end_index]/rc_shear_mean_history[0], 'r--', label='mean (w/o evidence)', linewidth=1.5 )
    ax2.plot(service_time[:end_index], rc_flex_cov_history[:end_index], 'b--', label='std (w/o evidence)', linewidth=1.5 )
    ax2.plot(service_time[:end_index], rc_shear_cov_history[:end_index], 'r--', label='std (w/o evidence)', linewidth=1.5 )
    # evidene 2 results
    ax1.plot(service_time[:end_index], rc_flex_mean_condition_state[:end_index]/rc_flex_mean_condition_state[0], 'b-', label='mean (w/o evidence)', linewidth=1.5 )
    ax1.plot(service_time[:end_index], rc_shear_mean_condition_state[:end_index]/rc_shear_mean_condition_state[0], 'r-', label='mean (w/o evidence)', linewidth=1.5 )
    ax2.plot(service_time[:end_index], rc_flex_cov_condition_state[:end_index], 'b-', label='std (w/o evidence)', linewidth=1.5 )
    ax2.plot(service_time[:end_index], rc_shear_cov_condition_state[:end_index], 'r-', label='std (w/o evidence)', linewidth=1.5 )
    # figure settings
    ax2.set_xlabel(r'service time (yr)')
    f.text(0.04, 0.5, r'Normalized resistance $R(t)/R_0$', ha='center', va='center', rotation='vertical')
    #ax1.grid(True)
    #ax2.grid(True)
    #ax1.legend(loc='lower left', prop={'size':12})
    #ax2.legend(loc='upper left', prop={'size':12})
    
    # arrow legend
    ax1.annotate('mean values', xy=(5, 0.727), bbox=dict(boxstyle="round", fc="w"))
    ax1.annotate('flexure (w/o evidence)', xy=(30.0, 0.988), xytext=(4, 0.864),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    ax1.annotate('flexure (condition state)', xy=(46, 0.962), xytext=(16.74, 0.810),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    ax1.annotate('shear (w/o evidence)', xy=(61, 0.922), xytext=(33.8, 0.77),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    ax1.annotate('shear (condition state)', xy=(72.2, 0.787), xytext=(63.1, 0.85),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    
    ax2.annotate('coefficient of variation', xy=(5, 0.173), bbox=dict(boxstyle="round", fc="w"))
    ax2.annotate('flexure (w/o evidence)', xy=(94, 0.138), xytext=(68.75, 0.15),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=1))
    ax2.annotate('flexure (condition state)', xy=(61.3, 0.133), xytext=(35.2, 0.146),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=1))
    ax2.annotate('shear (condition state)', xy=(74, 0.1667), xytext=(52, 0.176),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=1))
    ax2.annotate('shear (w/o evidence)', xy=(11, 0.1515), xytext=(8, 0.164),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=1))
    
    all_figures.append(f)
    
    ## show figures
    plt.show()
    
    # save figures
    isSave = raw_input('Save figures? (y/n)')
    if isSave == 'y':
        plt.close('all')
        fig_no = 1
        for fig in all_figures:
            fig.savefig(FIG_PATH+'fig'+str(fig_no)+'.eps')
            fig_no += 1
        print 'figures saveed. Postprocessing ends.'
    elif isSave == 'n':
        print 'Postprocessing ends.'
    else:
        print 'Illegal input: figures not saved'
    
