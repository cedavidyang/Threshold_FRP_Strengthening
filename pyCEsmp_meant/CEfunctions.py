# CE functions
import os
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

from constants import *
from pyCEsmp import *

def combine2smps(pfpre, Spf2pre, pfmain, Spf2main):
    
    COVpfmain = np.sqrt(Spf2main) / pfmain
    COVpfpre = np.sqrt(Spf2pre) / pfpre
    w = Spf2main / (Spf2main + Spf2pre)
    # w = COVpfmain / (COVpfmain + COVpfpre)
    pf = (1.0-w) * pfmain + w * pfpre
    Spf2 = (1.0-w)**2 * Spf2main + w**2 * Spf2pre
    COVpf = np.sqrt(Spf2) / pf
    
    return pf, Spf2, COVpf
    
def printResults(time_array, result_dict):
    print('{time:<10} {beta:<10} {pf:<10} {COVpf:<10} {pfpre:<10} {COVpfpre:<10} {pfmain:<10} {COVpfmain:<10}'.format(time='# time',
                       beta='beta', pf='pf', COVpf='COVpf', pfpre='pfpre', COVpfpre='COVpfpre', pfmain='pfmain', COVpfmain='COVpfmain'))
    for ti in time_array:
        pf = result_dict['pf'][time_array==ti]
        Spf2 = result_dict['Spf2'][time_array==ti]
        pfpre = result_dict['pfpre'][time_array==ti]
        Spf2pre = result_dict['Spf2pre'][time_array==ti]
        pfmain = result_dict['pfmain'][time_array==ti]
        Spf2main = result_dict['Spf2main'][time_array==ti]
        COVpf = np.sqrt(Spf2) / pf
        COVpfpre = np.sqrt(Spf2pre) / pfpre
        COVpfmain = np.sqrt(Spf2main) / pfmain
        beta = -stats.norm.ppf(pf)
        print('{time:<10d} {beta:<10.2f} {pf:<10.2e} {COVpf:<10.4f} {pfpre:<10.2e} {COVpfpre:<10.4f} {pfmain:<10.2e} {COVpfmain:<10.4f}'.format(time=int(ti),
                           beta=beta[0], pf=pf[0], COVpf=COVpf[0], pfpre=pfpre[0], COVpfpre=COVpfpre[0], pfmain=pfmain[0], COVpfmain=COVpfmain[0]))
                           

def comparePfs():
    time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT)
    time_array = np.insert(time_array, 0, 1.)
    
    #flex_pf = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'flexure', 'results_parallel.npz'))['pf']
    #shear_pf = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'shear', 'results_parallel.npz'))['pf']
    #sys_corr_pf = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'system_correlated', 'results_parallel.npz'))['pf']
    #sys_indp_pf = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'system_indp', 'results_parallel.npz'))['pf']
    
    flex_pf = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'with_evidence', 'flexure', 'results_parallel.npz'))['pf']
    shear_pf = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'with_evidence', 'shear', 'results_parallel.npz'))['pf']
    sys_corr_pf = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'with_evidence', 'system_correlated', 'results_parallel.npz'))['pf']
    sys_indp_pf = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'with_evidence', 'system_indp', 'results_parallel.npz'))['pf']
    
    sys_indp_pf_2 = 1. - (1.- flex_pf) * (1. - shear_pf)
    
    plt.close('all')
    plt.rc('font', family='serif', size=12)
    
    plt.semilogy(time_array, flex_pf, 'bv-', label="flexural resistance")
    plt.semilogy(time_array, shear_pf, 'g^-', label="shear resistance")
    plt.semilogy(time_array, sys_corr_pf, 'rs-', label="system resistance (correlated loads)")
    plt.semilogy(time_array, sys_indp_pf, 'yo-', label="system resistance (independent loads)")
    plt.semilogy(time_array, sys_indp_pf_2, 'yo--', label="system resistance (from component reliability)")
    plt.axhline(y=stats.norm.cdf(-4.2),color='k',ls='dashed')
    
    ax = plt.gca()
    ax.annotate(r'$\beta_{T}=4.2$ ($p_{f}=1.33 \times 10^{-5}$)', xy=(13.5, 1.8e-5), xytext=(13.5, 1.8e-5))
    ax.annotate('Structural repair required', xy=(63.5, 1.271e-5), xytext=(13., 0.0002116), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    
    #plt.plot(time_array, flex_pf, 'bv-', label="flexural resistance")
    #plt.plot(time_array, shear_pf, 'g^-', label="shear resistance")
    #plt.plot(time_array, sys_corr_pf, 'rs-', label="system resistance (correlated loads)")
    #plt.plot(time_array, sys_indp_pf, 'yo-', label="system resistance (independent loads)")
    #plt.plot(time_array, sys_indp_pf_2, 'yo--', label="system resistance (from component reliability)")

    plt.legend(loc='lower right', fontsize=12)
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'Failure probability')
             
def optimalSmpFunc(r_array, smpObject, ti):
    
    if smpObject.getSmpNum() != 1:
        smpObject.setSmpNum(1)
        
    rm_mean = r_array[0].stats(moments='m')
    rm_std = np.sqrt(r_array[0].stats(moments='v'))
    
    rv_mean = r_array[1].stats(moments='m')
    rv_std = np.sqrt(r_array[1].stats(moments='v'))
    
    rm_array = np.linspace(rm_mean-5.*rm_std, rm_mean+3.*rm_std, 100)
    rv_array = np.linspace(rv_mean-5.*rv_std, rv_mean+3.*rv_std, 100)
    
    X, Y = np.meshgrid(rm_array, rv_array)
    
    Z = np.array([])
    Z2 = np.array([])
    for y in rv_array:
        for x in rm_array:
            smps = np.array([[x, y]])
            tmps = smpObject.priorPDF(smps) * (1. - smpObject.condAvailability(smps, ti))
            Z = np.append(Z, tmps)
            
            tmps, dummy = smpObject.score_samples(smps)
            Z2 = np.append(Z2, tmps)
            
    Z.shape = X.shape
    Z2.shape = X.shape
    
    np.savez('samplingFunc_yr'+str(int(ti))+'.npz', rm_grid = X, rv_grid=Y, optimal=Z, adapted=Z2)  

    
def compareSmpFunc(ti):
    res_path = os.path.join(os.path.abspath('./'), 'data', 'reliability', 'with_evidence', 'system_indp')
    file_name = 'samplingFunc_yr'+str(int(ti))+'.npz'
    data = np.load(os.path.join(res_path, 'postprocessing', file_name))

    X = data['rm_grid']
    Y = data['rv_grid']
    Z = data['optimal']
    Z2 = data['adapted']
    
    plt.figure()
    CS1 = plt.contour(X, Y, Z)
    plt.figure()
    CS2 = plt.contour(X, Y, Z2)
    #plt.clabel(CS, inline=1, fontsize=10)
    #plt.title('Simplest default with labels')

def postprocessSmp(ti):

    time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT)
    time_array = np.insert(time_array, 0, 1.)
    
    ntime = np.array(ti).size
    res_multi_path = os.path.join(os.path.abspath('./'), 'data', 'reliability', 'no_evidence', 'system_indp')
    res_buch_path = os.path.join(os.path.abspath('./'), 'data', 'reliability', 'no_evidence',  'system_indp_1peak')
    
    plt.close('all')
    plt.rc('font', family='serif', size=12)
    # figure 1: compare pf
    res_multi = np.load( os.path.join(res_multi_path, 'results_parallel.npz') )['pf']
    res_buch = np.load( os.path.join(res_buch_path, 'results_parallel.npz') )['pf']
    plt.plot(time_array, res_multi, 'bo-',
             time_array, res_buch, 'rs-')
    ax = plt.gca()
    ax.annotate('Gaussian mixture', xy=(80, 1.2e-6), xytext=(30, 1.8e-6), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    ax.annotate('unimodal Gaussian', xy=(60, 4.5e-7), xytext=(8, 10.5e-7), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'Failure probability')
    
    f, axarr = plt.subplots(ntime, 3, sharex='col', sharey='row')
    # subfigure 1: optimal pf
    for itime in range(ntime):
        file_name = 'samplingFunc_yr'+str(int(ti[itime]))+'.npz'
        data = np.load(os.path.join(res_multi_path, 'postprocessing', file_name))
        data1 = data
        data2 = np.load(os.path.join(res_buch_path, 'postprocessing', file_name))
        X = data['rm_grid']
        Y = data['rv_grid']
        Z = data['optimal']
        Z1 = data1['adapted']
        Z2 = data2['adapted']
    
        axarr[itime,0].contour(X, Y, Z)
        axarr[itime,1].contour(X, Y, Z1)
        axarr[itime,1].annotate(r'$P_f$={:.2e}'.format(res_multi[time_array==ti[itime]][0]), xy=(1000, 1000), xytext=(1000, 1000))
        axarr[itime,2].contour(X, Y, Z2)
        axarr[itime,2].annotate(r'$P_f$={:.2e}'.format(res_buch[time_array==ti[itime]][0]), xy=(1000, 1000), xytext=(1000, 1000))
        if itime == 0:
            axarr[itime, 0].set_title( 'optimal',fontsize=12 )
            axarr[itime, 1].set_title( 'mixture',fontsize=12 )
            axarr[itime, 2].set_title( 'unimodal',fontsize=12 )
        f.text(0.95, 0.25+1./3.5*2 -1./3.5*(itime) , 'Year '+str(ti[itime]), ha='center', va='center', rotation='vertical')
    
    f.text(0.04, 0.5, 'Shear strength (kN)', ha='center', va='center', rotation='vertical')
    f.text(0.5, 0.03, 'Flexural strength (kNm)', ha='center', va='center')
    
    plt.show()
    
    
def plotPf():
    time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT)
    time_array = np.insert(time_array, 0, 1.)

    res_covt_path = os.path.join(os.path.abspath('./'), 'data', 'reliability', 'with_evidence', 'system_indp_covt')
    res_cov0_path = os.path.join(os.path.abspath('./'), 'data', 'reliability', 'with_evidence',  'system_indp_cov0')
    res_covmean_path = os.path.join(os.path.abspath('./'), 'data', 'reliability', 'with_evidence',  'system_indp_covmean')
    
    plt.close('all')
    plt.rc('font', family='serif', size=12)
    # figure 1: compare pf
    res_covt = np.load( os.path.join(res_covt_path, 'results_parallel.npz') )['pf']
    res_cov0 = np.load( os.path.join(res_cov0_path, 'results_parallel.npz') )['pf']
    res_covmean = np.load( os.path.join(res_covmean_path, 'results_parallel.npz') )['pf']
    plt.semilogy(time_array, res_covt, 'bo-',
             time_array, res_cov0, 'rs-',
             time_array, res_covmean, 'g^-')
    ax = plt.gca()
    #ax.annotate('Gaussian mixture', xy=(80, 1.2e-6), xytext=(30, 1.8e-6), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    #ax.annotate('unimodal Gaussian', xy=(60, 4.5e-7), xytext=(8, 10.5e-7), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'Failure probability')
    
    plt.show()    
    
    
def plotResistance():
    FLEXTURE_NO_EV_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data', 'no_evidence')
    SHEAR_NO_EV_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data', 'no_evidence_shear')
    FLEXTURE_WITH_EV_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_condition_state')
    SHEAR_WITH_EV_DATAFILE_PATH = os.path.join(os.path.abspath('./'), 'data', 'evidence_condition_state_shear')
    
    # read degradation data
    flexure_no_ev_datafile = os.path.join(FLEXTURE_NO_EV_DATAFILE_PATH, 'LWS_results.txt')
    shear_no_ev_datafile = os.path.join(SHEAR_NO_EV_DATAFILE_PATH, 'LWS_results.txt')
    flexure_with_ev_datafile = os.path.join(FLEXTURE_WITH_EV_DATAFILE_PATH, 'LWS_results.txt')
    shear_with_ev_datafile = os.path.join(SHEAR_WITH_EV_DATAFILE_PATH, 'LWS_results.txt')
    service_time = np.loadtxt(flexure_no_ev_datafile)[0,:]
    # flexural strength rc
    flexure_no_ev_mean_history = np.loadtxt(flexure_no_ev_datafile)[19,:]
    flexure_no_ev_std_history = np.loadtxt(flexure_no_ev_datafile)[20,:]
    flexure_with_ev_mean_history = np.loadtxt(flexure_with_ev_datafile)[19,:]
    flexure_with_ev_std_history = np.loadtxt(flexure_with_ev_datafile)[20,:]
    # shear strength rc
    shear_no_ev_mean_history = np.loadtxt(shear_no_ev_datafile)[21,:]
    shear_no_ev_std_history = np.loadtxt(shear_no_ev_datafile)[22,:]
    shear_with_ev_mean_history = np.loadtxt(shear_with_ev_datafile)[21,:]
    shear_with_ev_std_history = np.loadtxt(shear_with_ev_datafile)[22,:]
    
    def gd_func(x, a, b, c, d, e):
        return 1.0 + a*x +b*x**2+c*x**3+d*x**4+e*x**5
    gd_flex_popt_no_ev, pcov = curve_fit(gd_func, service_time, flexure_no_ev_mean_history/flexure_no_ev_mean_history[0])
    gd_shear_popt_no_ev, pcov = curve_fit(gd_func, service_time, shear_no_ev_mean_history/shear_no_ev_mean_history[0])
    gd_flex_popt_with_ev, pcov = curve_fit(gd_func, service_time, flexure_with_ev_mean_history/flexure_with_ev_mean_history[0])
    gd_shear_popt_with_ev, pcov = curve_fit(gd_func, service_time, shear_with_ev_mean_history/shear_with_ev_mean_history[0])
    
    gd_flexure_no_ev = lambda x: gd_func(x, gd_flex_popt_no_ev[0], gd_flex_popt_no_ev[1], gd_flex_popt_no_ev[2], gd_flex_popt_no_ev[3], gd_flex_popt_no_ev[4])
    gd_shear_no_ev = lambda x: gd_func(x, gd_shear_popt_no_ev[0], gd_shear_popt_no_ev[1], gd_shear_popt_no_ev[2], gd_shear_popt_no_ev[3], gd_shear_popt_no_ev[4])

    gd_flexure_with_ev = lambda x: gd_func(x, gd_flex_popt_with_ev[0], gd_flex_popt_with_ev[1], gd_flex_popt_with_ev[2], gd_flex_popt_with_ev[3], gd_flex_popt_with_ev[4])
    gd_shear_with_ev = lambda x: gd_func(x, gd_shear_popt_with_ev[0], gd_shear_popt_with_ev[1], gd_shear_popt_with_ev[2], gd_shear_popt_with_ev[3], gd_shear_popt_with_ev[4])

    # compare
    plt.close('all')
    plt.rc('font', family='serif', size=12)
    
    plt.plot(service_time, flexure_no_ev_mean_history/flexure_no_ev_mean_history[0], 'b.')
    plt.plot(service_time, gd_flexure_no_ev(service_time), 'b--')
    plt.plot(service_time, shear_no_ev_mean_history/shear_no_ev_mean_history[0], 'r.')
    plt.plot(service_time, gd_shear_no_ev(service_time), 'r--')
             
    plt.plot(service_time, flexure_with_ev_mean_history/flexure_with_ev_mean_history[0], 'b.')
    plt.plot(service_time, gd_flexure_with_ev(service_time), 'b-')
    plt.plot(service_time, shear_with_ev_mean_history/shear_with_ev_mean_history[0], 'r.')
    plt.plot(service_time, gd_shear_with_ev(service_time), 'r-')
    
    ax = plt.gca()
    ax.annotate('Updated shear resistance\n(data and polynomial fitting)', xy=(64, 0.775), xytext=(8, 0.682), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    ax.annotate('Prior shear', xy=(92, 0.89), xytext=(60, 0.81), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    ax.annotate('Updated flexure', xy=(82, 0.925), xytext=(44, 0.85), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    ax.annotate('Prior flexure', xy=(68, 0.97), xytext=(35, 0.895), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'Mean resistance / initial mean resistance')
    
    #plt.plot(service_time, flexure_no_ev_mean_history, 'bo')
    #plt.plot(service_time, gd_flexure_no_ev(service_time)*flexure_no_ev_mean_history[0], 'b--')
    #plt.plot(service_time, shear_no_ev_mean_history, 'rs')
    #plt.plot(service_time, gd_shear_no_ev(service_time)*shear_no_ev_mean_history[0], 'r--')
    #         
    #plt.plot(service_time, flexure_with_ev_mean_history, 'b^')
    #plt.plot(service_time, gd_flexure_with_ev(service_time)*flexure_with_ev_mean_history[0], 'b-')
    #plt.plot(service_time, shear_with_ev_mean_history, 'rv')
    #plt.plot(service_time, gd_shear_with_ev(service_time)*shear_with_ev_mean_history[0], 'r-')
    
    plt.show()
    
    
def postprocessStrengthening():
    time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT)
    time_array = np.insert(time_array, 0, 1.)
    
    sys_indp_pf = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'with_evidence', 'system_indp', 'results_parallel.npz'))['pf']
    pf_interp1d = interp1d(time_array, sys_indp_pf)
    
    sys_indp_pf_str1 = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'frp_U_anchor', 'with_evidence', 'system_indp', 'results_parallel.npz'))['pf']
    sys_indp_pf_str2 = np.load(os.path.join(os.path.abspath('./'), 'data', 'reliability', 'frp_U_bond', 'with_evidence', 'system_indp', 'results_parallel.npz'))['pf']
    
    tmp_time_array = np.arange(RELIABILITY_DT,SERVICE_LIFE+RELIABILITY_DT,RELIABILITY_DT) - EM_SHEAR_YR
    new_time_array = np.insert(tmp_time_array[tmp_time_array>0], 0, 1.)+EM_SHEAR_YR
    new_time_array = np.insert(new_time_array, 0, EM_SHEAR_YR)
    sys_indp_pf_str1 = np.hstack((pf_interp1d(EM_SHEAR_YR),sys_indp_pf_str1))
    sys_indp_pf_str2 = np.hstack((pf_interp1d(EM_SHEAR_YR),sys_indp_pf_str2))
    
    #indx = time_array>EM_SHEAR_YR
    #sys_indp_pf_part1 = sys_indp_pf[np.logical_not(indx)]
    #sys_indp_pf_part1 = np.append(sys_indp_pf_part1, pf_interp1d(EM_SHEAR_YR))
    #sys_indp_pf_part2 = 1-(1-sys_indp_pf[indx]) / (1-pf_interp1d(EM_SHEAR_YR))
    #sys_indp_pf_part2 = np.insert(sys_indp_pf_part2, 0, 1-(1-pf_interp1d(EM_SHEAR_YR+1)) / (1-pf_interp1d(EM_SHEAR_YR)))
    #sys_indp_pf = np.hstack((sys_indp_pf_part1,sys_indp_pf_part2))
    #time_array = np.hstack((time_array, EM_SHEAR_YR, EM_SHEAR_YR+1))
    #time_array.sort()
    
    plt.close('all')
    plt.rc('font', family='serif', size=12)
    
    plt.semilogy(time_array[time_array>=50], sys_indp_pf[time_array>=50], 'bo-', label="w/o strengthening")
    plt.semilogy(new_time_array[new_time_array>=50], sys_indp_pf_str1[new_time_array>=50], 'r^-', label="U-jacketing w/ anchor")
    plt.semilogy(new_time_array[new_time_array>=50], sys_indp_pf_str2[new_time_array>=50], 'gv-', label="U-jacketing w/o anchor")
    plt.axhline(y=stats.norm.cdf(-4.2),color='k',ls='dashed')

    ax = plt.gca()
    ax.annotate(r'$\beta_{T}=4.2$ ($p_{f}=1.33 \times 10^{-5}$)', xy=(52, 1.8e-5), xytext=(52, 1.8e-5))
    #ax.annotate('Structural repair required', xy=(63.5, 1.271e-5), xytext=(13., 0.0002116), arrowprops=dict(width=0.2, headwidth=8, facecolor='k', shrink=0.05))

    #plt.legend(loc='lower right', fontsize=12)
    plt.xlabel(r'service time (yr)')
    plt.ylabel(r'Failure probability')
    plt.xlim((48, 132))
    plt.ylim((1e-7, 1e-3))
 
    plt.show()
    
#if __name__ == '__main__':
    #postprocessStrengthening()
    #SERVICE_LIFE = 100
    #postprocessSmp(np.array([80, 90, 100]))
