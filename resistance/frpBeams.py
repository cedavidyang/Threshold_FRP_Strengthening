# FRP-stengthened beam samples

import numpy as np
import math
import scipy.optimize as opt
import scipy.special as spec
import scipy.stats as stats
from scipy.optimize import fsolve

from pyre import *
from resistance.beams import RCBeam

class FRPBeam(RCBeam):
    def __init__(self, name, nsmp, ME_dict, material_dict, geo_dict):
        """FRP-strengthened beam model
            :Attributs:
                h   --- beam height (mm)
                a   --- distance between FRP plate end and support (mm)
                ss  --- shear span (mm)
                l   --- span (mm)
                dc  --- compressive edge to centroid of compressive reinforcement (mm)
                Asc --- area of compressive reinforcement (mm2)
                ft  --- concrete tensile strength (MPa)
                fcuk --- characteristic value of concrete cube (MPa)
                ecu --- ultimate strain of concrete
                Es  --- modulus of steel (MPa)
                fyc --- yield strength of compressive reinforcement (mm2)
                Efrp --- modulus of FRP (MPa)
                ffrp --- strength of FRP (MPa)
                tfrp --- total thickness of FRP (mm)
                bfrp --- width of FRP (mm)
                fanchor --- contribution of anchor to FRP stress (MPa)
                eini --- initial strain under sevice load
                ME_flex  --- model error for rc flexural strength (aci model)
        """
        self.nsmp = nsmp
        # geometric data
        self.h = np.array(geo_dict['h'])
        self.a = np.array(geo_dict['a'])    
        self.ss = np.array(geo_dict['ss'])
        self.l = np.array(geo_dict['l'])
        self.dc = np.array(geo_dict['dc'])
        self.Asc =  np.array(geo_dict['Asc'])
        # material data
        self.ft = np.array(material_dict['ft'])
        self.fcuk = material_dict['fcuk']
        self.ecu = np.array(material_dict['ecu'])
        self.Es = np.array(material_dict['Es'])
        self.fyc = np.array(material_dict['fyc'])
        self.Efrp = np.array(material_dict['Efrp'])
        self.ffrp = np.array(material_dict['ffrp'])
        self.tfrp = np.array(material_dict['tfrp'])
        self.bfrp = np.array(material_dict['bfrp'])
        self.fanchor = np.array(material_dict['fanchor'])
        self.eini = np.array(material_dict['eini'])
        self.env = np.array(material_dict['env'])
        
        RCBeam.__init__(self, name, ME_dict, material_dict, geo_dict)
        
        # deduced parameters
        self.fcu = self.fc / 0.8
        self.Ec = 3460*np.sqrt(self.fcu) + 3210
        # 0.67 = 0.8 * 0.85 = cylinder/cube * in-situ/cylinder. if cylinder strength is used, 0.67 should be replaced by 0.85
        self.eco = 2*0.67*self.fcu / self.Ec 
        self.Afrp = self.tfrp * self.bfrp
        
        # for shear strengthening
        self.dfrp  = np.array(geo_dict['dfrp'])
        self.dfrpt  = np.array(geo_dict['dfrpt'])
        self.barType  = np.array(geo_dict['barType'])
        self.dsv  = np.array(geo_dict['dsv'])
        self.wfrpv  = np.array(geo_dict['wfrpv'])
        self.sfrpv  = np.array(geo_dict['sfrpv'])
        self.frpIncline = np.array(geo_dict['frpIncline'])
        
        self.Efrpv  = np.array(material_dict['Efrpv'])
        self.tfrpv  = np.array(material_dict['tfrpv'])
        self.ffrpv  = np.array(material_dict['ffrpv'])
        self.strform = material_dict['strform']
        
        
    def frpStress(self, x, sfmax, ismp):
        """get stress in FRP flexural strengthening
           :Args:
               x --- location of neutral axis
               sfmax --- maximum stress in FRP, either controlled by FRP rupture or IC debonding
        """

        # determine FRP stress
        e_fe = self.ecu[ismp] * (self.h[ismp] - x) / x - self.eini[ismp]
        sfrp = self.Efrp[ismp] * e_fe
        sfrp = np.minimum(sfrp, sfmax)
        
        return sfrp
        
        
    def tensionSteelStress(self, x, sfrp, ismp):
        efrp = sfrp / self.Efrp[ismp]
        e_s = (efrp+self.eini[ismp]) * (self.d[ismp]-x) / (self.h[ismp]-x)
        tmp = self.Es[ismp] * e_s
        ssteel = np.minimum(tmp, self.fy[ismp])
        
        return ssteel
        
        
    def compSteelStress(self, x, sfrp, ismp):
        efrp = sfrp / self.Efrp[ismp]
        e_sc = (x-self.dc[ismp]) / (self.h[ismp]-x) * (efrp+self.eini[ismp])
        ssteel = self.Es[ismp] * e_sc
        if ssteel>self.fyc[ismp]:
            ssteel = self.fyc[ismp]
        elif ssteel<-self.fyc[ismp]:
            ssteel = -self.fyc[ismp]
        
        return ssteel
        
        
    def forceBalance(self, x, sfmax, ismp):
        """ get force balance residual
           :Args:
               x --- location of neutral axis
               sfmax --- maximum stress in FRP, either controlled by FRP rupture or IC debonding
        """
        
        sfrp = self.frpStress(x, sfmax, ismp)
        fs= self.tensionSteelStress(x, sfrp, ismp)
        fsc = self.compSteelStress(x, sfrp, ismp)
        
        # calculate concrete strain for FRP control cases. From subfunction FRP_stress, the following equation equals to e_cu for concrete control cases.
        efrp = sfrp / self.Efrp[ismp]
        e_cf = (efrp + self.eini[ismp] ) * x / (self.h[ismp]-x)
        if e_cf>0 and e_cf<=self.eco[ismp]:
            k1 = ( self.Ec[ismp] * e_cf / (2*0.8*self.fcu[ismp])) * (1.- self.Ec[ismp]*e_cf / (6*0.8*self.fcu[ismp]))
        elif e_cf>self.eco[ismp] and e_cf<=self.ecu[ismp]:
            k1 = 1. - self.eco[ismp] / 3. / e_cf
        
        if self.bf[ismp] == self.b[ismp]:    # rectangular beam
            Ac = self.b[ismp] * x
        elif x<=self.hf[ismp] and self.bf[ismp] != self.b[ismp]:    # type I T-beam
            Ac = self.bf[ismp] * x
        elif x>self.hf[ismp] and self.bf[ismp] != self.b[ismp]:    # type II T-beam
            Ac = self.b[ismp] * x + (self.bf[ismp]-self.b[ismp]) * self.hf[ismp]
        r = self.As[ismp]*fs + self.Afrp[ismp]*sfrp - self.Asc[ismp]*fsc - k1*(0.8*self.fcu[ismp])*Ac

        return r
        
        
    def flexCapacity(self):
        """flexural capacity of FRP-strengthened beams (HK model), neglect initial strain
        """        
        e_ini = self.eini
        
        # determine debonding strain of FRP       
        beta_w = np.sqrt( (2.-self.bfrp/self.b) / (1.+self.bfrp/self.b) )
        # beta_w = np.sqrt( (2.25-bfrp/b) / (1.25+bfrp/b) )
        t_max = 1.5*beta_w*self.ft
        Lee = 0.228*np.sqrt(self.Efrp*self.tfrp)
        Ld = self.ss - self.a
        alpha = 3.41*Lee/Ld
        f_anch = self.fanchor 
        fdbic = 0.114*(4.41-alpha)*t_max*np.sqrt(self.Efrp/self.tfrp) + f_anch
        ff = np.minimum(self.ffrp, fdbic)
        
        # determine location of neutral axis
        Mu = np.zeros(self.nsmp)
        for ismp in range(self.nsmp):
            c_ini = 0.1*self.d(ismp)
            func = lambda x: self.forceBalance(x, ff[ismp], ismp)
            c, infodict, ier, msg = fsolve( func, c_ini)
            if ier != 1:
                print 'convergence failed'

            sfrp = self.frpStress(c, ff[ismp])
            fs = self.tensionSteelStress(c, sfrp)
            fsc = self.compSteelStress(c, sfrp)
            
            e_fe = sfrp / self.Efrp[ismp]
            e_ini = self.eini[ismp]
            e_cf = (e_fe + e_ini ) * c / ( self.h[ismp]-c)
            if e_cf > self.ecu[ismp]:
                e_cf = self.ecu[ismp]
            k2 = 0.33+0.045*e_cf/self.eco[ismp]
            Mu[ismp] = 1e-6 * ( self.As[ismp] * fs * (self.d[ismp]-k2*c)
                              + self.Afrp[ismp] * sfrp * (self.h[ismp]-k2*c)
                              + self.Asc[ismp] * fsc * (k2*c-self.dc[ismp]) )
            Mu = self.ME_flex * Mu
        
        return Mu
        
        
    def shearCapacity(self):
        # maximum frp stress
        zt = self.dfrpt
        zb = (self.dv-(self.h-self.dfrp))-0.1*self.dv
        h_frp_e = zb-zt
        
        beta = self.frpIncline
        ratio = self.wfrpv / (self.sfrpv * np.sin(beta))
        temp = (2.-ratio) / (1.+ratio)
        temp[temp<0] = 0.
        beta_w = np.sqrt(temp)
        
        if self.strform.lower() == 'w' or self.strform.lower() =='u':
            l_max = h_frp_e / np.sin(beta)
        elif self.strform.lower() == 'side':
            l_max = h_frp_e / 2. / np.sin(beta)

        l_e = np.sqrt( self.Efrpv * self.tfrpv / np.sqrt(0.8*self.fcu) )
        lmbd = l_max / l_e
        
        beta_l = np.zeros(self.nsmp)
        beta_l[lmbd >= 1] = 1.0
        beta_l[lmbd < 1 ] = np.sin(np.pi*lmbd[lmbd<1]/2.)
        
        sdb = self.env * 0.315 * beta_w * beta_l * np.sqrt(self.Efrpv * np.sqrt(0.8*self.fcu) / self.tfrpv)
        
        failMode = np.zeros(self.nsmp, dtype=int)
        if self.strform.lower() == 'w':
            s_frp_max = 0.8*self.ffrpv
            failMode = np.ones(self.nsmp)    # 1 for fiber rupture
        else:
            s_frp_max = np.minimum(sdb, 0.8*self.ffrpv)
            failMode[sdb > 0.8*self.ffrpv] = 1
            failMode[sdb <= 0.8*self.ffrpv] = 2    # 2 for debonding failure
    
        # distribution factor and FRP contribution
        D_frp = np.zeros(self.nsmp)
        zeta = zt/zb
        D_frp[failMode==1] = (1. + zeta[failMode==1]) / 2.
        indx = np.logical_and(failMode==2, lmbd>1)
        D_frp[indx] = 1. - (np.pi-2.) / (np.pi*lmbd[indx])
        indx = np.logical_and(failMode==2, lmbd<=1)
        D_frp[indx] = 2. / np.pi / lmbd[indx] * (1.-np.cos(np.pi/2.*lmbd[indx])) / np.sin(np.pi/2.*lmbd[indx])
        
        f_frp_e = s_frp_max * D_frp
        V_frp = 2. * f_frp_e * self.tfrpv * self.wfrpv * h_frp_e * (np.sin(beta)+np.cos(beta)) / self.sfrpv / 1e3
        V_frp[V_frp < 0] = 0.
        
        # steel contribution and steel-FRP interaction
        if self.strform.lower() == 'side':
            phi_s = np.zeros(self.nsmp)
            A = np.zeros(self.nsmp)
            Vs = np.zeros(self.nsmp)
            mu = np.zeros(self.nsmp)
            kfrp = np.zeros(self.nsmp)
            
            if self.barType.lower() == 'plain':
                phi_s = 1e5 / (self.dsv**1.13 * self.fyv**1.71)
                A = phi_s * (2.045*2*np.sin(beta)*lmbd+3.24)
            elif self.barType.lower() == 'deformed':
                phi_s = 1e5 / (self.dsv**0.834 * self.fyv**1.88)
                A = phi_s * (1.01*2*np.sin(beta)*lmbd+2.13)
            
            Vs = 1e-3*2*np.pi*self.dsv**2/4. * self.fyv * self.dv / self.sv
            Vs[Vs<0] = 0
            
            mu = Vs / V_frp
            mu[V_frp == 0] = 0
            kfrp = A / (A + mu)
            
        else:
            Vs = 1e-3*2*np.pi*self.dsv**2/4. * self.fyv * self.dv / self.sv
            Vs[Vs<0] = 0
            kfrp = np.ones(self.nsmp)
            
                
        V_frp_interaction = 2*kfrp * f_frp_e * self.tfrpv * self.wfrpv * h_frp_e * (np.sin(beta)+np.cos(beta)) / self.sfrpv / 1e3
        V_frp_interaction[V_frp_interaction<0] = 0
        
        # concrete contribution
        temp1 = (400. / self.dv)**0.25
        
        if self.fcuk<=40:
            vr = 0.4*np.ones(self.nsmp)
        elif self.fcuk>40 and self.fcuk<80:
            vr = 0.4*(self.fcuk/40.)**(2./3.)
        elif self.fcuk>=80:
            vr = 0.4*(80./40.)**(2./3.)
            
        As_min = vr * self.b * self.sv / (0.87*self.fyv)
        As = self.Asvt
        temp1[temp1<0.8] = 0.8 
        temp1[ np.logical_and(As>=As_min, temp1<1) ] = 1.0
        
        Vc = 0.79*(100*0.01)**(1./3.) * temp1 * self.b * self.dv * 1e-3
        
        if self.fcuk>25 and self.fcuk<80:
            Vc = Vc * (self.fcuk/25)**(1./3.)
        elif self.fcuk>=80:
            Vc = Vc * (80./25.)**(1./3.)
            
        Vc[Vc<0] = 0.0
        
        Vtotal = V_frp_interaction + Vs + Vc
        temp3 = np.minimum(np.sqrt(self.fcu), 7*1.25/(1-1.645*FC_COV))
        indx = Vtotal*1e3 / (self.b*self.dv) > temp3
        Vtotal[indx]= temp3[indx] * self.b[indx] * self.dv[indx] / 1e3
        Vtotal = self.ME_shear * Vtotal
        
        return Vtotal
        

if __name__ == '__main__':
    from constants import *
    from resistance.resistanceFuncs import *
    from corrosion import *
    
#    ## check 1, Monte Carlo simulation
#    seed_indx = np.arange(1,55,2)
#    ## seeds for sample generation
#    CONCRETE_STRENGTH_SEED = seed_indx[6]**2
#    CONCRETE_MODULUS_SEED = seed_indx[7]**2
#    CONCRETE_TENSION_SEED = seed_indx[8]**2
#    # flexural resistance
#    ME_FLEX_RC_SEED = seed_indx[12]**2
#    FSY_SEED = seed_indx[13]**2
#    BEAM_WIDTH_SEED = seed_indx[14]**2
#    BEAM_DEPTH_SEED = seed_indx[15]**2
#    FLANGE_WIDTH_SEED = seed_indx[16]**2
#    FLANGE_DEPTH_SEED = seed_indx[17]**2
#    # shear resistance
#    ME_SHEAR_RC_SEED = seed_indx[18]**2
#    SHEAR_DEPTH_SEED = seed_indx[19]**2
#    FSYV_SEED = seed_indx[20]**2
#    SHEAR_INTERVAL_SEED = seed_indx[21]**2
#    # frp-related
#    ME_FLEX_FRP_SEED = seed_indx[22]**2
#    ME_SHEAR_FRP_SEED = seed_indx[23]**2
#    EFRP_SEED = seed_indx[23]**2
#    FFRP_SEED = seed_indx[24]**2
#    EFRPV_SEED = seed_indx[25]**2
#    FFRPV_SEED = seed_indx[26]**2
#
#    ## initial samples
#
#    # concrete strength
#    compressive_strength = concStrengthVariable()
#    np.random.seed(CONCRETE_STRENGTH_SEED)
#    fc_smp = compressive_strength.rv.rvs(size = N_SMP)
#        
#    # effective concrete elastic modulus
#    elastic_modulus = concEffEcVariable()
#    np.random.seed(CONCRETE_MODULUS_SEED)
#    Ec_smp = elastic_modulus.rv.rvs(size = N_SMP)
#        
#    # concrete tensile strength
#    tensile_strength = concTensileVariable()
#    np.random.seed(CONCRETE_TENSION_SEED)
#    ft_smp = tensile_strength.rv.rvs(size = N_SMP)
#        
#    # model error of flexural strength of RC beams
#    ME_flex = modelErrorRCFlexVariable()
#    np.random.seed(ME_FLEX_RC_SEED)
#    ME_flex_smp = ME_flex.rv.rvs(size=N_SMP)
#        
#    # yielding strength of steel
#    fy_steel = steelYieldingVariable()
#    np.random.seed(FSY_SEED)
#    fy_smp = fy_steel.rv.rvs(size=N_SMP)
#        
#    # beam width and depth
#    beam_width = beamWidthVariable()
#    np.random.seed(BEAM_WIDTH_SEED)
#    b_smp = beam_width.rv.rvs(size=N_SMP)
#    beam_depth = beamDepthVariable()
#    np.random.seed(BEAM_DEPTH_SEED)
#    d_smp = beam_depth.rv.rvs(size=N_SMP)
#        
#    # flange width and depth
#    flange_width = flangeWidthVariable()
#    np.random.seed(FLANGE_WIDTH_SEED)
#    bf_smp = flange_width.rv.rvs(size=N_SMP)
#    flange_depth = flangeDepthVariable()
#    np.random.seed(FLANGE_DEPTH_SEED)
#    hf_smp = flange_depth.rv.rvs(size=N_SMP)
#        
#    # model error of shear strength of RC beams
#    ME_shear = modelErrorRCShearVariable()
#    np.random.seed(ME_SHEAR_RC_SEED)
#    ME_shear_smp = ME_shear.rv.rvs(size=N_SMP)
#        
#    # yielding strength of shear reinforcement
#    fyv_steel = shearYieldingVariable()
#    np.random.seed(FSYV_SEED)
#    fyv_smp = fyv_steel.rv.rvs(size=N_SMP)
#        
#    # shear depth
#    dv = shearDepthVariable()
#    np.random.seed(SHEAR_DEPTH_SEED)
#    dv_smp = dv.rv.rvs(size=N_SMP)
#        
#    # shear interval
#    sv = shearIntervalVariable()
#    np.random.seed(SHEAR_INTERVAL_SEED)
#    sv_smp = sv.rv.rvs(size=N_SMP)       
#        
#    ## FRP-related random variables
#    # model error of flexural strengthening
#    ME_flex_frp = modelErrorFRPFlexVariable()
#    np.random.seed(ME_FLEX_FRP_SEED)
#    ME_flex_frp_smp = ME_flex_frp.rv.rvs(size=N_SMP)
#    # model error of shear strengthening
#    ME_shear_frp = modelErrorFRPShearVariable()
#    np.random.seed(ME_SHEAR_FRP_SEED)
#    ME_shear_frp_smp = ME_shear_frp.rv.rvs(size=N_SMP)
#    # material properties
#    ecu_smp = ECU_MEAN * np.ones(N_SMP)
#    Es_smp = ES_MEAN * np.ones(N_SMP)
#    fyc_smp = FYC_MEAN * np.ones(N_SMP)
#    # flexural strengthening
#    Efrp = EfrpVariable()
#    np.random.seed(EFRP_SEED)
#    Efrp_smp = Efrp.rv.rvs(size=N_SMP)
#        
#    ffrp = ffrpVariable()
#    np.random.seed(FFRP_SEED)
#    ffrp_smp = ffrp.rv.rvs(size=N_SMP)
#        
#    tfrp_smp = TFRP_MEAN * np.ones(N_SMP)
#    bfrp_smp = BFRP_MEAN * np.ones(N_SMP)
#        
#    if EM_FLEX_FORM.lower() == 'ic':
#        fanchor_smp = np.zeros(N_SMP)
#    else:
#        fanchor_smp = FANCHOR_MEAN * np.ones(N_SMP)
#    
#    eini_smp = np.zeros(N_SMP)
#        
#    # shear strengthening
#    Efrpv = EfrpShearVariable()
#    np.random.seed(EFRPV_SEED)
#    Efrpv_smp = Efrpv.rv.rvs(size=N_SMP)
#    
#    ffrpv = ffrpShearVariable()
#    np.random.seed(FFRPV_SEED)
#    ffrpv_smp = ffrpv.rv.rvs(size=N_SMP)
#    
#    tfrpv_smp = TFRPV_MEAN * np.ones(N_SMP) 
#    
#    # geometric properties
#    h_smp = d_smp +  BEAM_HEIGHT_MEAN - BEAM_DEPTH_MEAN
#    a_smp = PLATE2SUPPORT_MEAN * np.ones(N_SMP)
#    ss_smp = SHEAR_SPAN_MEAN * np.ones(N_SMP)
#    l_smp = SPAN_MEAN * np.ones(N_SMP)
#    dsc_smp = DSC_MEAN * np.ones(N_SMP)
#    Asc_smp = np.pi/4*dsc_smp**2 * DSC_NO
#    dfrp_smp = np.copy(h_smp)
#    dfrpt_smp = hf_smp
#    wfrpv_smp = WFRPV_MEAN * np.ones(N_SMP)
#    sfrpv_smp = SFRPV_MEAN * np.ones(N_SMP)
#    frpIncline_smp = FRP_INCLINE_MEAN * np.ones(N_SMP) * np.pi/180.
#
#    
#    ds_smp = np.ones(N_SMP) * DS_REGION1_MEAN
#    ds_smp2 = np.ones(N_SMP) * DS_REGION2_MEAN
#    ds_prev = np.ones(N_SMP) * DS_REGION1_MEAN
#    ds_prev2 = np.ones(N_SMP) * DS_REGION2_MEAN
#
#
#    # get resistance
#    # name, ME_flex, fc, fy, Ast, b, d, bf, hf
#    # flexural
#    Ast_smp = np.pi/4 * ds_smp**2 * (REGION1_NO + REGION2_NO)
#    ME_dict, material_dict, geo_dict = assembleBeamDict(ME_flex_smp, ME_shear_smp,
#                                                        fc_smp, fy_smp, fyv_smp,
#                                                        Ast_smp, Ast_smp, b_smp, d_smp, bf_smp, hf_smp, dv_smp, sv_smp)
#    rcBeam = RCBeam('rc_beam', ME_dict, material_dict, geo_dict)
#    mu_smp = rcBeam.flexCapacity()
#    # shear
#    Asvt_smp = np.pi/4 * ds_smp**2 * (REGION1_NO + REGION2_NO)
#    ME_dict, material_dict, geo_dict = assembleBeamDict(ME_flex_smp, ME_shear_smp,
#                                                        fc_smp, fy_smp, fyv_smp,
#                                                        Asvt_smp, Asvt_smp, b_smp, d_smp, bf_smp, hf_smp, dv_smp, sv_smp)
#                                                        
#    rcBeam = RCBeam('rc_beam', ME_dict, material_dict, geo_dict)
#    vu_smp = rcBeam.shearCapacity()   
#        
#    ME_dict, material_dict, geo_dict = assembleFrpBeamDict(ME_dict, material_dict, geo_dict, ME_flex_frp_smp, ME_shear_frp_smp,
#            ft_smp, FC_NOM/0.8, ecu_smp, Es_smp, fyc_smp, Efrp_smp, ffrp_smp, tfrp_smp, bfrp_smp, fanchor_smp, eini_smp,
#            Efrpv_smp, tfrpv_smp, ffrpv_smp, EM_SHEAR_FORM,
#            h_smp, a_smp, ss_smp, l_smp, dsc_smp, Asc_smp, dfrp_smp, dfrpt_smp, BAR_TYPE, ds_smp, wfrpv_smp, sfrpv_smp, frpIncline_smp)
#    frpBeam = FRPBeam('frp_beam', N_SMP, ME_dict, material_dict, geo_dict)
#    vufrp_smp = frpBeam.shearCapacity()

    ## compare with hand calculation
    nsmp = 1
    # concrete strength
    fc_smp = np.array([FC_NOM])  
    ## concrete elastic modulus
    #Ec_smp = 4733.0 * np.sqrt(FC_NOM)    
    # concrete tensile strength
    ft_smp = np.array([0.6227 * np.sqrt(FC_NOM)])
    # model error of flexural strength of RC beams
    ME_flex_smp = np.array([1.0])
    # yielding strength of steel
    fy_smp = np.array([FSY_NOM])      
    # beam width and depth
    b_smp = np.array([WEB_WIDTH_NOM])
    d_smp = np.array([BEAM_DEPTH_NOM])
    # flange width and depth
    bf_smp = np.array([FLANGE_WIDTH_NOM])
    hf_smp = np.array([FLANGE_THICK_NOM])
    # model error of shear strength of RC beams
    ME_shear_smp = np.array([1.0])
    # yielding strength of shear reinforcement
    fyv_smp = np.array([FSYV_NOM])
    # shear depth
    dv_smp = np.array([BEAM_DEPTH_NOM])
    # shear interval
    sv_smp = np.array([SHEAR_INTERVAL_NOM])     
    ## FRP-related random variables
    # model error of flexural strengthening
    ME_flex_frp_smp = np.array([1.0])
    # model error of shear strengthening
    ME_shear_frp_smp = np.array([1.0])
    # material properties
    ecu_smp = np.array([ECU_NOM])
    Es_smp = np.array([ES_NOM])
    fyc_smp = np.array([FYC_NOM])
    # flexural strengthening
    Efrp_smp = np.array([EFRP_NOM])
    ffrp_smp = np.array([FFRP_NOM])
    tfrp_smp = np.array([TFRP_MEAN])
    bfrp_smp = np.array([BFRP_MEAN])
        
    if EM_FLEX_FORM.lower() == 'ic':
        fanchor_smp = np.array([0.0])
    else:
        fanchor_smp = np.array([FANCHOR_NOM])
    
    eini_smp = np.array([0.0])
    
    # shear strengthening
    Efrpv_smp = np.array([EFRPV_NOM])
    ffrpv_smp = np.array([FFRPV_NOM])
    tfrpv_smp = np.array([TFRPV_NOM])
    # geometric properties
    h_smp = d_smp +  BEAM_HEIGHT_NOM - BEAM_DEPTH_NOM
    a_smp = np.array([PLATE2SUPPORT_NOM])
    ss_smp = np.array([SHEAR_SPAN_NOM])
    l_smp = np.array([SPAN_NOM])
    dsc_smp = np.array([DSC_NOM])
    Asc_smp = np.pi/4*dsc_smp**2 * DSC_NO
    dfrp_smp = h_smp
    dfrpt_smp = hf_smp
    wfrpv_smp = np.array([WFRPV_NOM])
    sfrpv_smp = np.array([SFRPV_NOM])
    frpIncline_smp = np.array([FRP_INCLINE_NOM]) * np.pi/180.
    ds_smp = np.array([DS_REGION1_NOM])
    ds_smp2 = np.array([DS_REGION2_NOM])
    ds_prev = np.array([DS_REGION1_NOM])
    ds_prev2 = np.array([DS_REGION2_NOM])

    # get resistance
    # shear
    Asvt_smp = np.pi/4 * ds_smp**2 * (REGION1_NO + REGION2_NO)
    ME_dict, material_dict, geo_dict = assembleBeamDict(ME_flex_smp, ME_shear_smp,
                                                        fc_smp, fy_smp, fyv_smp,
                                                        Asvt_smp, Asvt_smp, b_smp, d_smp, bf_smp, hf_smp, dv_smp, sv_smp) 
    ME_dict, material_dict, geo_dict = assembleFrpBeamDict(ME_dict, material_dict, geo_dict, ME_flex_frp_smp, ME_shear_frp_smp,
            ft_smp, FC_NOM/0.8, ecu_smp, Es_smp, fyc_smp, Efrp_smp, ffrp_smp, tfrp_smp, bfrp_smp, fanchor_smp, eini_smp,
            Efrpv_smp, tfrpv_smp, ffrpv_smp, EM_SHEAR_FORM,
            h_smp, a_smp, ss_smp, l_smp, dsc_smp, Asc_smp, dfrp_smp, dfrpt_smp, BAR_TYPE, ds_smp*10.4/12.7, wfrpv_smp, sfrpv_smp, frpIncline_smp)
    frpBeam = FRPBeam('frp_beam', nsmp, ME_dict, material_dict, geo_dict)
    vufrp_smp = frpBeam.shearCapacity()