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
    datafile = os.path.join(FIG_DATAFILE_PATH, 'no_evidence', 'LWS_results.txt')    
    service_time = np.loadtxt(datafile)[0,:]
    
    yrOfInterest = 100.0
    end_index = np.where(service_time==yrOfInterest)[0]+1
      
    # chloride history
    chloride_mean_history = np.loadtxt(datafile)[1,:]
    chloride_std_history = np.loadtxt(datafile)[2,:]
    
    # corrosion initiation prob
    corrosion_prob_history = np.loadtxt(datafile)[3,:]
    
    # corrosion rate history
    corrosion_rate_mean_history = np.loadtxt(datafile)[4,:]
    corrosion_rate_std_history = np.loadtxt(datafile)[5,:]
    
    # mean corrosion rate history (outdated no use)
    mean_corrosion_rate_mean_history = np.loadtxt(datafile)[6,:]
    mean_corrosion_rate_std_history = np.loadtxt(datafile)[7,:]
    
    # residual diameter history
    residual_diameter_mean_history = np.loadtxt(datafile)[8,:]
    residual_diameter_std_history = np.loadtxt(datafile)[9,:]
    
    # radial pressure history
    radial_pressure_mean_history = np.loadtxt(datafile)[10,:]
    radial_pressure_std_history = np.loadtxt(datafile)[11,:]
    
    # cracking initiation prob
    crack_initiation_history = np.loadtxt(datafile)[12,:]
    
    # ds at crack history
    ds_crack_mean_history = np.loadtxt(datafile)[13,:]
    ds_crack_std_history = np.loadtxt(datafile)[14,:]
    
    # crack width history
    crack_width_mean_history = np.loadtxt(datafile)[15,:]
    crack_width_std_history = np.loadtxt(datafile)[16,:]
    
    # diffusion coef. increase (fw)
    crack_diffusion_mean_history = np.loadtxt(datafile)[17,:]
    crack_diffusion_std_history = np.loadtxt(datafile)[18,:]
    
    
    ## =======================================================================##
    ## get results with evidence
    ## =======================================================================##
    datafile = os.path.join(FIG_DATAFILE_PATH, 'evidence_ini_crk_shear', 'LWS_results.txt')
      
    # chloride LWS
    chloride_mean_LWS = np.loadtxt(datafile)[1,:]
    chloride_std_LWS = np.loadtxt(datafile)[2,:]
    
    # corrosion initiation prob
    corrosion_prob_LWS = np.loadtxt(datafile)[3,:]
    
    # corrosion rate LWS
    corrosion_rate_mean_LWS = np.loadtxt(datafile)[4,:]
    corrosion_rate_std_LWS = np.loadtxt(datafile)[5,:]
    
    # mean corrosion rate LWS (outdated no use)
    mean_corrosion_rate_mean_LWS = np.loadtxt(datafile)[6,:]
    mean_corrosion_rate_std_LWS = np.loadtxt(datafile)[7,:]
    
    # residual diameter LWS
    residual_diameter_mean_LWS = np.loadtxt(datafile)[8,:]
    residual_diameter_std_LWS = np.loadtxt(datafile)[9,:]
    
    # radial pressure LWS
    radial_pressure_mean_LWS = np.loadtxt(datafile)[10,:]
    radial_pressure_std_LWS = np.loadtxt(datafile)[11,:]
    
    # cracking initiation prob
    crack_initiation_LWS = np.loadtxt(datafile)[12,:]
    
    # ds at crack LWS
    ds_crack_mean_LWS = np.loadtxt(datafile)[13,:]
    ds_crack_std_LWS = np.loadtxt(datafile)[14,:]
    
    # crack width LWS
    crack_width_mean_LWS = np.loadtxt(datafile)[15,:]
    crack_width_std_LWS = np.loadtxt(datafile)[16,:]
    
    # diffusion coef. increase (fw)
    crack_diffusion_mean_LWS = np.loadtxt(datafile)[17,:]
    crack_diffusion_std_LWS = np.loadtxt(datafile)[18,:]
    

    ## =======================================================================##
    ## generate and save figures
    ## =======================================================================##
    
    plt.close('all')
    plt.rc('font', family='serif', size=12)
    all_figures = []
    
    # results  of Cl
    f = plt.figure()
    plt.hold(True)
    # MC reults
    plt.plot(service_time[:end_index], chloride_mean_history[:end_index], 'b-', label='mean (w/o evidence)', linewidth=0.5 )
    plt.plot(service_time[:end_index], chloride_std_history[:end_index], 'b--', label='std (w/o evidence)', linewidth=0.5 )
    # LWS results
    plt.plot(service_time[:end_index], chloride_mean_LWS[:end_index], 'r-', label='mean (w/ evidence)', linewidth=2)
    plt.plot(service_time[:end_index], chloride_std_LWS[:end_index], 'r--', label='std (w/ evidence)', linewidth=2)
    # figure settings
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'chloride content ($kg/m^3$)')
    plt.grid(True)
    plt.legend(loc='upper left', prop={'size':12})
    all_figures.append(f)
    
    
    # icorr
    f = plt.figure()
    plt.hold(True)
    # MC reults
    #plt.plot(service_time[:end_index], corrosion_rate_mean_history[:end_index], 'b-', label='mean (w/o evidence)', linewidth=1 )
    #plt.plot(service_time[:end_index], corrosion_rate_std_history[:end_index], 'b--', label='std (w/o evidence)', linewidth=1 )
    # LWS results
    plt.plot(service_time[:end_index], corrosion_rate_mean_LWS[:end_index], 'b-', label='mean (w/ evidence)' , linewidth=1.5)
    plt.plot(service_time[:end_index], corrosion_rate_std_LWS[:end_index], 'b--', label='std (w/ evidence)' , linewidth=1.5)
    plt.axvline(x=20., color='k', ls='--')
    plt.axvline(x=52., color='k', ls='--')
    ax = plt.gca()
    ax.annotate('Phase I', xy = (4, 0.45), bbox=dict(boxstyle="round", fc="w"))
    ax.annotate(r'Phase II-III', xy = (24, 0.45), bbox=dict(boxstyle="round", fc="w"))
    ax.annotate(r'Phase IV-', xy = (56, 0.45), bbox=dict(boxstyle="round", fc="w"))
    ax.annotate('corrosion initiation', xy=(20, 0.18), xytext=(25, 0.23), arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    ax.annotate('crack initiation', xy=(52, 0.18), xytext=(57, 0.23), arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    ax.annotate(r'mean $i_{corr}$', xy=(60, 0.31), xytext=(70, 0.37), arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    ax.annotate(r'std $i_{corr}$', xy=(60, 0.11), xytext=(70, 0.06), arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3', facecolor='k', shrinkB=0))
    # figure settings
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'corrosion rate $i_{corr}$ ($\mu A/cm^2$)')
    plt.grid(False)
    #plt.legend(loc='upper right', prop={'size':12})
    all_figures.append(f)
    
    
    '''# imean (outdated , no use)
    f = plt.figure()
    plt.hold(True)
    # MC reults
    plt.plot(service_time[:end_index], mean_corrosion_rate_mean_history[:end_index], 'b-', label='mean (w/o evidence)' , linewidth=0.5)
    plt.plot(service_time[:end_index], mean_corrosion_rate_std_history[:end_index], 'b--', label='mean (w/o evidence)' , linewidth=0.5)
    # LWS results
    plt.plot(service_time[:end_index], mean_corrosion_rate_mean_LWS[:end_index], 'r-', label='mean (w/ evidence)' , linewidth=2)
    plt.plot(service_time[:end_index], mean_corrosion_rate_std_LWS[:end_index], 'r--', label='std (w/ evidence)' , linewidth=2)
    # figure settings
    plt.rc('font', family='serif', size=12)
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'corrosion rate $i_{corr}$ ($\mu A/cm^2$)')
    plt.grid(True)
    plt.legend(loc='lower right', prop={'size':12})
    all_figures.append(f)'''
    
    # residual diameter mean
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    # MC reults
    ax1.plot(service_time[:end_index], residual_diameter_mean_history[:end_index], 'b-', label='mean (w/o evidence)', linewidth=0.5 )
    ax2.plot(service_time[:end_index], residual_diameter_std_history[:end_index], 'b--', label='std (w/o evidence)', linewidth=0.5 )
    ## LWS results
    ax1.plot(service_time[:end_index], residual_diameter_mean_LWS[:end_index], 'r-', label='mean (w/ evidence)' , linewidth=2)
    ax2.plot(service_time[:end_index], residual_diameter_std_LWS[:end_index], 'r--', label='std (w/ evidence)' , linewidth=2)
    # figure settings
    ax2.set_xlabel(r'service time (yr)')
    f.text(0.04, 0.5, r'residual diameter (mm)', ha='center', va='center', rotation='vertical')
    ax1.grid(True)
    ax2.grid(True)
    ax1.legend(loc='lower left', prop={'size':12})
    ax2.legend(loc='lower right', prop={'size':12})
    all_figures.append(f)
    
    # radial pressure
    f = plt.figure()
    plt.hold(True)
    # MC reults
    plt.plot(service_time[:end_index], radial_pressure_mean_history[:end_index], 'b-', label='mean (w/o evidence)', linewidth=0.5 )
    plt.plot(service_time[:end_index], radial_pressure_std_history[:end_index], 'b--', label='std (w/o evidence)', linewidth=0.5 )
    # LWS results
    plt.plot(service_time[:end_index], radial_pressure_mean_LWS[:end_index], 'r-', label='mean (w/ evidence)' , linewidth=2)
    plt.plot(service_time[:end_index], radial_pressure_std_LWS[:end_index], 'r--', label='std (w/ evidence)' , linewidth=2)
    # figure settings
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'radial pressure (MPa)')
    plt.grid(True)
    plt.legend(loc='upper left', prop={'size':12})
    all_figures.append(f)
    
    # plot the figure
    f = plt.figure()
    plt.hold(True)
    # MC reults
    plt.plot(service_time[:end_index], crack_initiation_history[:end_index], 'b-', label='crack initiation (MC)', linewidth=0.5)
    plt.plot(service_time[:end_index], corrosion_prob_history[:end_index], 'b--', label='corrosion initiation (MC)', linewidth=0.5)
    # LWS results
    plt.plot(service_time[:end_index], crack_initiation_LWS[:end_index], 'r-', label='crack initiation (LWS)', linewidth = 2)
    plt.plot(service_time[:end_index], corrosion_prob_LWS[:end_index], 'r--', label='corrosion initiation (LWS)', linewidth = 2)
    # figure settings
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'CDF')
    plt.grid(True)
    plt.legend(loc='upper left', prop={'size':12})
    all_figures.append(f)
    
    '''# plot the figure
    f = plt.figure()
    plt.hold(True)
    # MC reults
    plt.plot(service_time[:end_index], ds_crack_mean_history[:end_index], 'b-', label='mean (w/o evidence)', linewidth=0.5 )
    plt.plot(service_time[:end_index], ds_crack_std_history[:end_index], 'b--', label='std (w/o evidence)', linewidth=0.5 )
    # LWS results
    plt.plot(service_time[:end_index], ds_crack_mean_LWS[:end_index], 'r-', label='mean (w/ evidence)' , linewidth=2)
    plt.plot(service_time[:end_index], ds_crack_std_LWS[:end_index], 'r--', label='std (w/ evidence)' , linewidth=2)
    # figure settings
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'ds at crack (mm)')
    plt.grid(True)
    plt.legend(loc='upper left', prop={'size':12})
    all_figures.append(f)'''
    
    
    # crack width
    f = plt.figure()
    plt.hold(True)
    # MC reults
    plt.plot(service_time[:end_index], crack_width_mean_history[:end_index], 'b-', label='mean (w/o evidence)', linewidth=0.5 )
    plt.plot(service_time[:end_index], crack_width_std_history[:end_index], 'b--', label='std (w/o evidence)', linewidth=0.5 )
    # LWS results
    plt.plot(service_time[:end_index], crack_width_mean_LWS[:end_index], 'r-', label='mean (w/ evidence)' , linewidth=2)
    plt.plot(service_time[:end_index], crack_width_std_LWS[:end_index], 'r--', label='std (w/ evidence)' , linewidth=2)
    # figure settings
    plt.rc('font', family='serif', size=12)
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'crack width (mm)')
    plt.grid(True)
    plt.legend(loc='upper left', prop={'size':12})
    all_figures.append(f)
    
    # diffusion of cracked concrete
    f = plt.figure()
    plt.hold(True)
    # MC reults
    plt.plot(service_time[:end_index], crack_diffusion_mean_history[:end_index], 'b-', label='mean (w/o evidence)', linewidth=0.5 )
    plt.plot(service_time[:end_index], crack_diffusion_std_history[:end_index], 'b--', label='std (w/o evidence)', linewidth=0.5 )
    # LWS results
    plt.plot(service_time[:end_index], crack_diffusion_mean_LWS[:end_index], 'r-', label='mean (w/ evidence)' , linewidth=2)
    plt.plot(service_time[:end_index], crack_diffusion_std_LWS[:end_index], 'r--', label='std (w/ evidence)' , linewidth=2)
    # figure settings
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'diffusion of cracked concrete ($D_{cracked} / D_{uncracked}$)')
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
            fig.savefig(FIG_PATH+'fig'+str(fig_no)+'.eps')
            fig_no += 1
        print 'figures saveed. Postprocessing ends.'
    elif isSave == 'n':
        print 'Postprocessing ends.'
    else:
        print 'Illegal input: figures not saved'
    