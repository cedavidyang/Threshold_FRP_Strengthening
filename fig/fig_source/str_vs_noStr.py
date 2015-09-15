# check cl and cdf of corrosion initiation cdf
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

def main():
    # FIG_PATH = os.path.join(os.path.abspath('./'), 'fig')
    FIG_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data')
    dataNoStr = os.path.join(FIG_DATAFILE_PATH,
                    'evidence_condition_state_shear', 'LWS_results.txt')
    dataStr = os.path.join(FIG_DATAFILE_PATH, 'frp_U_anchor_degrade',
                  'evidence_condition_state_shear', 'LWS_results.txt')
    service_time = np.loadtxt(dataNoStr)[0,:]

    yrOfInterest = 130.0
    end_index = np.where(service_time==yrOfInterest)[0]+1

    # residual shear resistance without strengthening
    rc_shear_mean_history = np.loadtxt(dataNoStr)[21,:]
    rc_shear_std_history = np.loadtxt(dataNoStr)[22,:]
    rc_shear_cov_history = rc_shear_std_history/rc_shear_mean_history
    # residual shear resistance with strengthening
    frp_shear_mean_history = np.loadtxt(dataStr)[21,:]
    frp_shear_std_history = np.loadtxt(dataStr)[22,:]
    frp_shear_cov_history = frp_shear_std_history/frp_shear_mean_history

    # plot figures
    plt.close('all')
    plt.rc('font', family='serif', size=12)

    f, ax1 = plt.subplots(1)
    ax1.plot(service_time[:end_index],
             rc_shear_mean_history[:end_index]/rc_shear_mean_history[0],
             'b-', label='mean (w/o strengthening)', linewidth=1.5 )
    ax1.plot(service_time[:end_index],
             frp_shear_mean_history[:end_index]/frp_shear_mean_history[0],
             'r-', label='mean (w/ strengthening)', linewidth=1.5 )
    plt.show()

if __name__ == '__main__':
    main()
