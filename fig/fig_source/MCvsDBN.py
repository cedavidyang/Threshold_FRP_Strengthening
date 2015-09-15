# check cl and cdf of corrosion initiation cdf
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

FIG_PATH = os.path.join(os.path.abspath('./'), 'fig')  
FIG_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data')   
    
if __name__ == '__main__':
    ## =======================================================================##
    ## get MC results
    ## =======================================================================##
    datafile = os.path.join(FIG_DATAFILE_PATH, 'LWS_results.txt')
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
    ## get DBN results
    ## =======================================================================##
    datafile = FIG_DATAFILE_PATH+'DBN_results.txt' 
    
    # chloride DBN    
    chloride_mean_DBN = np.loadtxt(datafile)[0, :]
    chloride_std_DBN = np.loadtxt(datafile)[1, :]
    
    # corrosion initiation prob DBN
    corrosion_prob_DBN = np.loadtxt(datafile)[2,:]
    
    # corrosion rate DBN
    corrosion_rate_mean_DBN = np.loadtxt(datafile)[4,:]
    corrosion_rate_std_DBN = np.loadtxt(datafile)[5,:]
    
    # mean corrosion rate DBN (outdated no use)
    mean_corrosion_rate_mean_DBN = np.loadtxt(datafile)[6,:]
    mean_corrosion_rate_std_DBN = np.loadtxt(datafile)[7,:]
    
    # residual diameter DBN
    residual_diameter_mean_DBN = np.loadtxt(datafile)[8,:]
    residual_diameter_std_DBN = np.loadtxt(datafile)[9,:]
    
    # radial pressure DBN
    radial_pressure_mean_DBN = np.loadtxt(datafile)[10,:]
    radial_pressure_std_DBN = np.loadtxt(datafile)[11,:]
    
    # crack initiation probs
    crack_initiation_DBN = np.loadtxt(datafile)[12,:]
    
    # ds at cracking DBN
    ds_crack_mean_DBN = np.loadtxt(datafile)[14, :]
    ds_crack_std_DBN = np.loadtxt(datafile)[15, :]
    
    # crack width DBN
    crack_width_mean_DBN = np.loadtxt(datafile)[16,:]
    crack_width_std_DBN = np.loadtxt(datafile)[17,:]
    
    # diffusion coef. increase (fw)
    crack_diffusion_mean_DBN10 = np.loadtxt(datafile)[18,:]
    crack_diffusion_std_DBN10 = np.loadtxt(datafile)[19,:]
    

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
    plt.plot(service_time[:end_index], chloride_mean_history[:end_index], 'b-', label='mean (MC)' )
    plt.plot(service_time[:end_index], chloride_std_history[:end_index], 'b--', label='std (MC)' )
    # DBN results
    plt.plot(service_time[:end_index], chloride_mean_DBN[:end_index], 'r-', label='mean (DBN)' , linewidth=2)
    plt.plot(service_time[:end_index], chloride_std_DBN[:end_index], 'r--', label='std (DBN)' , linewidth=2)
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
    plt.plot(service_time[:end_index], corrosion_rate_mean_history[:end_index], 'b-', label='mean (MC)' )
    plt.plot(service_time[:end_index], corrosion_rate_std_history[:end_index], 'b--', label='std (MC)' )
    # DBN results
    plt.plot(service_time[:end_index], corrosion_rate_mean_DBN[:end_index], 'r-', label='mean (DBN)' , linewidth=2)
    plt.plot(service_time[:end_index], corrosion_rate_std_DBN[:end_index], 'r--', label='std (DBN)' , linewidth=2)
    # figure settings
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'corrosion rate $i_{corr}$ ($\mu A/cm^2$)')
    plt.grid(True)
    plt.legend(loc='lower right', prop={'size':12})
    all_figures.append(f)
    
    
    '''# imean (outdated , no use)
    f = plt.figure()
    plt.hold(True)
    # MC reults
    plt.plot(service_time[:end_index], mean_corrosion_rate_mean_history[:end_index], 'b-', label='mean (MC)' )
    plt.plot(service_time[:end_index], mean_corrosion_rate_std_history[:end_index], 'b--', label='mean (MC)' )
    # DBN results
    plt.plot(service_time[:end_index], mean_corrosion_rate_mean_DBN[:end_index], 'r-', label='mean (DBN)' , linewidth=2)
    plt.plot(service_time[:end_index], mean_corrosion_rate_std_DBN[:end_index], 'r--', label='std (DBN)' , linewidth=2)
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
    ax1.plot(service_time[:end_index], residual_diameter_mean_history[:end_index], 'b-', label='mean (MC)' )
    ax2.plot(service_time[:end_index], residual_diameter_std_history[:end_index], 'b--', label='std (MC)' )
    ## DBN results
    #ax1.plot(service_time[:end_index], residual_diameter_mean_DBN[:end_index], 'r-', label='mean (DBN)' , linewidth=2)
    #ax2.plot(service_time[:end_index], residual_diameter_std_DBN[:end_index], 'r--', label='std (DBN)' , linewidth=2)
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
    plt.plot(service_time[:end_index], radial_pressure_mean_history[:end_index], 'b-', label='mean (MC)' )
    plt.plot(service_time[:end_index], radial_pressure_std_history[:end_index], 'b--', label='std (MC)' )
    # DBN results
    plt.plot(service_time[:end_index], radial_pressure_mean_DBN[:end_index], 'r-', label='mean (DBN)' , linewidth=2)
    plt.plot(service_time[:end_index], radial_pressure_std_DBN[:end_index], 'r--', label='std (DBN)' , linewidth=2)
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
    plt.plot(service_time[:end_index], crack_initiation_history[:end_index], 'b-', label='crack initiation (MC)')
    plt.plot(service_time[:end_index], corrosion_prob_history[:end_index], 'b--', label='corrosion initiation (MC)')
    # DBN results
    plt.plot(service_time[:end_index], crack_initiation_DBN[:end_index], 'r-', label='crack initiation (DBN)', linewidth = 2)
    plt.plot(service_time[:end_index], corrosion_prob_DBN[:end_index], 'r--', label='corrosion initiation (DBN)', linewidth = 2)
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
    plt.plot(service_time[:end_index], ds_crack_mean_history[:end_index], 'b-', label='mean (MC)' )
    plt.plot(service_time[:end_index], ds_crack_std_history[:end_index], 'b--', label='std (MC)' )
    # DBN results
    plt.plot(service_time[:end_index], ds_crack_mean_DBN[:end_index], 'r-', label='mean (DBN)' , linewidth=2)
    plt.plot(service_time[:end_index], ds_crack_std_DBN[:end_index], 'r--', label='std (DBN)' , linewidth=2)
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
    plt.plot(service_time[:end_index], crack_width_mean_history[:end_index], 'b-', label='mean (MC)' )
    plt.plot(service_time[:end_index], crack_width_std_history[:end_index], 'b--', label='std (MC)' )
    # DBN results
    plt.plot(service_time[:end_index], crack_width_mean_DBN[:end_index], 'r-', label='mean (DBN)' , linewidth=2)
    plt.plot(service_time[:end_index], crack_width_std_DBN[:end_index], 'r--', label='std (DBN)' , linewidth=2)
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
    plt.plot(service_time[:end_index], crack_diffusion_mean_history[:end_index], 'b-', label='mean (MC)' )
    plt.plot(service_time[:end_index], crack_diffusion_std_history[:end_index], 'b--', label='std (MC)' )
    # DBN results
    plt.plot(service_time[:end_index], crack_diffusion_mean_DBN10[:end_index], 'r-', label='mean (DBN)' , linewidth=2)
    plt.plot(service_time[:end_index], crack_diffusion_std_DBN10[:end_index], 'r--', label='std (DBN)' , linewidth=2)
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
    