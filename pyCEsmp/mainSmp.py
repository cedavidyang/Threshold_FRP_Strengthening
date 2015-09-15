# cross-entropy-based adaptive importance sampling
import numpy as np
import scipy.stats as stats

from scipy.integrate import quad
from sklearn.mixture.gmm import GMM
from constants import *
from pyre import *

import warnings

class MainSmp(object):
    def __init__(self, nsmp=None, gmdistr=None, rv_array=None, r_array=None, sl_array=None, gd_array=None, gcov_array=None, preSmp=None):
        """Create object:
        :Args:
            nsmp: number of samples
            gmdistr: Gaussian mixture distribution object
            r_array: numpy array of resistance [dtype = distr object]
            sl_array: numpy array of live load [dtype = double], assumed to be deterministic
            gd_array: numpy array of degradation functions [dtype = func object]
        """
        if preSmp==None:
            if any((nsmp==None, gmdistr==None, rv_array==None, r_array==None, sl_array==None, gd_array==None, gcov_array==None)):
                print '[warning: not enough parameters provided to construct MainSmp]'
            else:
                self.nsmp = nsmp
                self.gmdistr = gmdistr
                self.rv_array = rv_array
                self.r_array = r_array
                self.sl_array = sl_array
                self.gd_array = gd_array
                self.gcov_array = gcov_array
        
        else:
            self.nsmp = preSmp.nsmp
            self.gmdistr = preSmp.gmdistr
            self.rv_array = preSmp.rv_array
            self.r_array = preSmp.r_array
            self.sl_array = preSmp.sl_array
            self.gd_array = preSmp.gd_array
            self.gcov_array = preSmp.gcov_array
  
#=============================================================================#  
    def setSmpNum(self, new_nsmp):
        self.nsmp = new_nsmp
        
    def getSmpNum(self):
        return self.nsmp
        
    def setGmdistr(self, new_gmdistr):
        self.gmdistr = new_gmdistr
        
    def getGmdistr(self):
        return self.gmdistr
        
    def setRarray(self, new_r, indx=None):
        if indx == None:
            self.r_array = new_r
        else:
            self.r_array[indx] = new_r
        
    def getRarray(self):
        return self.r_array
        
    def setSLarray(self, new_sl, indx=None):
        if indx == None:
            self.sl_array = new_sl
        else:
            self.sl_array[indx] = new_sl
        
    def getSLarray(self):
        return self.sl_array
        
    def setGdarray(self, new_gd, indx=None):
        if indx == None:
            self.gd_array = new_gd
        else:
            self.gd_array[indx] = new_gd
        
    def getGdarray(self):
        return self.gd_array
#=============================================================================#
    
    def sample(self):
        smps = self.gmdistr.sample(n_samples=self.nsmp)
        return smps
        
    def score_samples(self, smps):
        logprobs, responsibilities = self.gmdistr.score_samples(smps)
        pdfs = np.exp(logprobs)
        return pdfs, responsibilities
        
    def priorPDF(self, smps):
        """ uncertainty of dead load is neglected. flexural and shear resistance are independent. Different girders are independent. """
        pdfs = np.zeros(smps.shape)
        for i in range(smps.shape[1]):
            pdfs[:,i] = self.rv_array[i].pdf(smps[:,i])
            
        jointPdf = np.prod(pdfs, axis=1)
        
        return jointPdf            
    
    def condAvailability(self, smps, ti):
        # if more than one girder, reshape smps first
        
        # one girder
        gm = self.gd_array[0]
        gv = self.gd_array[1]
        gcovm = self.gcov_array[0]
        gcovv = self.gcov_array[1]
        sl = self.sl_array[0]
        
        sdm = M_DCDW_MEAN
        sdv = V_DCDW_MEAN
        cm = 1.0*LDF_FLEX
        # cv = V_DCDW_MEAN/M_DCDW_MEAN*LDF_SHEAR
        cv = V_LLIM_MEAN/M_LLIM_MEAN*LDF_SHEAR
        
        available = np.array([])
        for i in range(self.nsmp):
            #rm = smps[i,0]
            #rv = smps[i,1]
            rm = self.r_array[0]
            rv = self.r_array[1]
            covm = np.exp(smps[i,0])
            covv = np.exp(smps[i,1])
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    #ccdf_smps = lambda t: self.failAtT(t, rm, rv)
                    #ccdf_smps = np.vectorize(ccdf_smps)
                    
                    # # serial system with correlated Ml and Vl
                    # #print "SystemWithCorrelated"
                    # ccdf_smps = lambda t: 1. - sl.cdf( np.minimum( (rm*gm(t)*(covm**gcovm(t))-sdm)/cm,  (rv*gv(t)*(covv**gcovv(t))-sdv)/cv ) )
                
                    # serial system with independent Ml and Vl
                    #print "SystemWithIndependent"
                    ccdf_smps = lambda t: 1. - sl.cdf( (rm*gm(t)*(covm**gcovm(t))-sdm)/cm ) * sl.cdf( (rv*gv(t)*(covv**gcovv(t))-sdv)/cv )
                
                    # # component Ml
                    # #print "flexure"
                    # ccdf_smps = lambda t: 1. - sl.cdf( (rm*gm(t)*(covm**gcovm(t))-sdm)/cm )
                    
                    ## component Vl
                    ##print "shear"
                    #ccdf_smps = lambda t: 1. - sl.cdf( (rv*gv(t)*(covv**gcovv(t))-sdv)/cv )
                
            
                    tmp = np.exp( -LL_ARRIVAL_RATE * quad(ccdf_smps, 0, ti)[0] )
                    available = np.append(available, tmp)
                except Warning:
                    print 'failed integration'
            # available.append( np.exp( -LL_ARRIVAL_RATE * quad(ccdf_smps, 0, ti)[0] ) )
        
        return np.array(available)
        
        
    def failAtT(self, t, rm, rv):
        # one girder
        gm = self.gd_array[0]
        gv = self.gd_array[1]
        sl = self.sl_array[0]
        
        sdm = M_DCDW_MEAN
        sdv = V_DCDW_MEAN
        cm = 1.0*LDF_FLEX
        cv = V_LLIM_MEAN/M_LLIM_MEAN*LDF_SHEAR
        
        slt_mean = sl.stats(moments='m') * (1.+LL_MEAN_INCREASE_RATE)**t
        slt_std = np.sqrt(sl.stats(moments='v')) * (1.+LL_STD_INCREASE_RATE)**t
        slt = stats.norm(loc=slt_mean, scale=slt_std)
        
        ## serial system with correlated Ml and Vl
        ##print "SystemWithCorrelated"
        #ccdf_smps = lambda t: 1. - slt.cdf( np.minimum( (rm*gm(t)-sdm)/cm,  (rv*gv(t)-sdv)/cv ) )
        
        # serial system with independent Ml and Vl
        #print "SystemWithIndependent"
        ccdf = 1. - slt.cdf( (rm*gm(t)-sdm)/cm ) * slt.cdf( (rv*gv(t)-sdv)/cv )
        
        ## component Ml
        ##print "flexure"
        #ccdf_smps = lambda t: 1. - slt.cdf( (rm*gm(t)-sdm)/cm )
        
        ## component Vl
        ##print "shear"
        #ccdf_smps = lambda t: 1. - slt.cdf( (rv*gv(t)-sdv)/cv )
        
        return ccdf
        
        
    def getPf(self,smps, ti):
        f = self.priorPDF(smps)
        g = 1 - self.condAvailability(smps, ti)
        hv, responsibilities = self.score_samples(smps)
        gf2hv = g*f/hv 
        
        pfmain = np.mean(gf2hv)
        Spf2main = 1./self.nsmp * np.var(gf2hv, ddof=1)
        
        return pfmain, Spf2main
        
         
if __name__ == '__main__':
    ti = 1
    nCompo = 4
    nmain = 200
    muVopt = np.array([[891.404727222421, 913.501588431349, 905.090608103064, 920.997600729747],
                       [1828.76493195385, 1709.76493998393, 1729.75265083924, 1689.58888405073]]).T
    piVopt = np.array([0.591259881095148, 0.0904722135483567, 0.310287166921472, 0.00798073843502256])
    covarVopt = np.array([[[3637.0, -892.0],[-892.0, 95600.0]], [[8332.0, -6257.0],[-6257.0, 52611.0]],
                         [[6567.0, -4470.0],[-4470.0, 56876.0]], [[9948.0, -7620.0], [-7620.0, 48147.0]]])

    gmdistr = GMM(n_components=4, covariance_type='full')
    gmdistr.weights_ = piVopt
    gmdistr.means_ = muVopt
    gmdistr.covars_ = covarVopt
    resistRv = Lognormal('r', mean=1790.1, stdv=283.37)
    rDistr = resistRv.rv
    #r_array = np.array([stats.norm(loc=1790.1, scale=283.37), stats.norm(loc=1790.1, scale=283.37)])
    r_array = np.array([rDistr, rDistr])
    sl_array = np.array([stats.norm(loc=301.1, scale=120.44)])
    ti_mean = 15.91
    a =  0.0135;    b = 0.8580
    gd = lambda x: (1-a*(x-ti_mean)**b) if x>ti_mean else 1.0
    gd_array = np.array([gd, gd])
    
    testSmp = MainSmp(nmain, gmdistr, r_array, sl_array, gd_array)
    smps = testSmp.sample()
    f = testSmp.priorPDF(smps)
    g = 1-testSmp.condAvailability(smps, ti)
    hv, responsibilities = testSmp.score_samples(smps)
    
    pf, Spf2 = testSmp.getPf(smps, ti)
