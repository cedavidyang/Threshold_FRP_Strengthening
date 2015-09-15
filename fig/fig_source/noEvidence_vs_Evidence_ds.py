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
    residual_diameter_mean_history = np.loadtxt(datafile)[8,:]
    residual_diameter_std_history = np.loadtxt(datafile)[9,:]

    
    ## =======================================================================##
    ## get results with evidence
    ## =======================================================================##
    #datafile = FIG_DATAFILE_PATH+'evidence_ini_crk/LWS_results.txt' 
    ## residual diameter LWS
    #residual_diameter_mean_ini_crk = np.loadtxt(datafile)[8,:]
    #residual_diameter_std_ini_crk = np.loadtxt(datafile)[9,:]
    
    #datafile = FIG_DATAFILE_PATH+'evidence_condition_state/LWS_results.txt' 
    datafile = os.path.join(FIG_DATAFILE_PATH, 'evidence_condition_state_shear', 'LWS_results.txt')
    # residual diameter LWS
    residual_diameter_mean_condition_state= np.loadtxt(datafile)[8,:]
    residual_diameter_std_condition_state = np.loadtxt(datafile)[9,:]
    
    #datafile = FIG_DATAFILE_PATH+'evidence_half_cell/LWS_results.txt' 
    datafile = os.path.join(FIG_DATAFILE_PATH, 'evidence_half_cell_shear', 'LWS_results.txt')
    # residual diameter LWS
    residual_diameter_mean_half_cell = np.loadtxt(datafile)[8,:]
    residual_diameter_std_half_cell = np.loadtxt(datafile)[9,:]
    
    #datafile = FIG_DATAFILE_PATH+'evidence_corrosion_rate/LWS_results.txt'
    datafile = os.path.join(FIG_DATAFILE_PATH, 'evidence_corrosion_rate_shear', 'LWS_results.txt') 
    # residual diameter LWS
    residual_diameter_mean_corrosion_rate = np.loadtxt(datafile)[8,:]
    residual_diameter_std_corrosion_rate = np.loadtxt(datafile)[9,:]

    ## =======================================================================##
    ## generate and save figures
    ## =======================================================================##
    
    plt.close('all')
    plt.rc('font', family='serif', size=12)
    all_figures = []
    
    # residual diameter mean
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    # MC reults
    ax1.plot(service_time[:end_index], residual_diameter_mean_history[:end_index], 'b-', label='mean (w/o evidence)', linewidth=1.5 )
    ax2.plot(service_time[:end_index], residual_diameter_std_history[:end_index], 'b-', label='std (w/o evidence)', linewidth=1.5 )
    ##evidene 1 results
    #ax1.plot(service_time[:end_index], residual_diameter_ini_crk[:end_index], 'r-', label='mean (w/ evidence)' , linewidth=2)
    #ax2.plot(service_time[:end_index], residual_diameter_ini_crk[:end_index], 'r--', label='std (w/ evidence)' , linewidth=2)
    # evidene 2 results
    ax1.plot(service_time[:end_index], residual_diameter_mean_condition_state[:end_index], 'r--', label='mean (condition state)' , linewidth=1.5)
    ax2.plot(service_time[:end_index], residual_diameter_std_condition_state[:end_index], 'r--', label='std (condition state)' , linewidth=1.5)
    ## evidene 3 results
    ax1.plot(service_time[:end_index], residual_diameter_mean_half_cell[:end_index], 'g-.', label='mean (half-cell)' , linewidth=1.5)
    ax2.plot(service_time[:end_index], residual_diameter_std_half_cell[:end_index], 'g-.', label='std (half-cell)' , linewidth=1.5)
    ## evidene 4 results
    ax1.plot(service_time[:end_index], residual_diameter_mean_corrosion_rate[:end_index], 'y:', label='mean (corrosion rate)' , linewidth=1.5)
    ax2.plot(service_time[:end_index], residual_diameter_std_corrosion_rate[:end_index], 'y:', label='std (corrosion rate)' , linewidth=1.5)
    # figure settings
    ax2.set_xlabel(r'service time (yr)')
    f.text(0.04, 0.5, r'residual diameter (mm)', ha='center', va='center', rotation='vertical')
    #ax1.grid(True)
    #ax2.grid(True)
    #ax1.legend(loc='lower left', prop={'size':12})
    #ax2.legend(loc='upper left', prop={'size':12})
    
    # arrow legend
    ax1.annotate('mean (w/o evidence)', xy=(51.5, 12.1), xytext=(8.87, 11.0), arrowprops=dict(width=0.2, headwidth=8, facecolor='k'))
    ax1.annotate('mean (condition state)', xy=(73, 10.5), xytext=(29.87, 9.6), arrowprops=dict(width=0.2, frac=0.12, headwidth=8, facecolor='k'))
    ax1.annotate('mean (half-cell)', xy=(36.8, 12.4), xytext=(2.6, 11.4), arrowprops=dict(width=0.2, headwidth=8, facecolor='k')) 
    ax1.annotate('mean (corrosion rate)', xy=(63.1, 11.5), xytext=(14.8, 10.2), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))    

    ax2.annotate('std (w/o evidence)', xy=(75.2, 0.58), xytext=(28.6, 0.7), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    ax2.annotate('std (condition state)', xy=(86.5, 0.8), xytext=(40, 0.9), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    ax2.annotate('std (half-cell)', xy=(50.4, 0.229), xytext=(7.6, 0.37), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05)) 
    ax2.annotate('std (corrosion rate)', xy=(69.76, 0.44), xytext=(18, 0.56), arrowprops=dict(width=0.2, frac=0.08, headwidth=8, facecolor='k', shrink=0.05))    
    
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
    