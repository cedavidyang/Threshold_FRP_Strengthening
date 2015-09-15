# cross-entropy-based adaptive importance sampling
import numpy as np
import scipy.stats as stats

from scipy.integrate import quad
from sklearn.mixture.gmm import GMM
from constants import *

class MainSmp(object):
    def __init__(self, nmain, gmdistr, r_array, sl_array, gd_array):
        """Create object:
        :Args:
            nmain: number of samples
            gmdistr: Gaussian mixture distribution object
            r_array: numpy array of resistance [dtype = distr object]
            sl_array: numpy array of live load [dtype = double], assumed to be deterministic
            gd_array: numpy array of degradation functions [dtype = func object]
        """
        self.nmain = nmain
        self.gmdistr = gmdistr
        self.r_array = r_array
        self.sl_array = gd_array
  
#=============================================================================#  
    def setSmpNum(self, new_nmain):
        self.nmain = new_nmain
        
    def getSmpNum(self):
        return self.nmain
        
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
        smps = self.gmdistr.sample(n_samples=nmain)
        return smps
        
    def score_sample(self, smps):
        logprobs, responsibilities = self.gmdistr.score_sample(smps)
        return logprobs, responsibilities
        
    def priorPDF(self, smps):
        """ uncertainty of dead load is neglected. flexural and shear resistance are independent. Different girders are independent. """
        pdfs = np.zeros(smps.shape)
        for i in range(smps.shape[1]):
            pdfs[:,i] = r_array[i].pdf(smps[:,i])
            
        jointPdf = np.prod(pdfs, axis=1)
        
        return jointPdf            
    
    def condAvailability(self, smps, ti):
        # if more than one girder, reshape smps first
        
        # one girder
        gm = self.gd_array[0]
        gv = self.gd_array[1]
        sl = self.sl_array[0]
        sdm = M_DCDW_MEAN
        sdv = V_DCDW_MEAN
        cm = 1.0
        cv = V_DCDW_MEAN/M_DCDW_MEAN
        
        available = []
        for i in range(self.nmain):
            rm = smps[i,0]
            rv = smps[i,1]

            ccdf_smps = lambda t: 1 - sl.cdf( np.minimum( (rm*gm(t)-sdm)/cm,  (rv*gv(t)-sdv)/cv ) )
            ccdf_smps = np.vectorize(ccdf_smps)
            available.append( np.exp( -l_arr * quad(ccdf_smps, 0, ti) ) )
        
        return np.array(available)