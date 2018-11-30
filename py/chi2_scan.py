import numpy as np
from scipy.interpolate import interp2d

basedir = '/Users/jfarr/Projects/H0LyA/data/BOSS_scans/'
scan_locations =   {'cf': basedir+'/BOSSDR11LyaF_k.scan',
                    'xcf': basedir+'/BOSSDR11QSOLyaF.scan'
                    }

class chi2_interpolators():
    def __init__(self,scan_locations,DA_over_rd_fid,c_over_Hrd_fid,kind='linear'):
        #Create a dictionary containing an interpolator for each scan.
        interpolators = {}
        for corr_type in scan_locations:
            scan = np.loadtxt(scan_locations[corr_type])
            interpolators[corr_type] = interp2d(scan[:,0],scan[:,1],scan[:,2],kind='linear')
        #Add the dictionary to the object.
        self.interpolators = interpolators
        self.DA_over_rd_fid = DA_over_rd_fid
        self.c_over_Hrd_fid = c_over_Hrd_fid
        return
    def get_chi2_distances(self,DA_over_rd,c_over_Hrd,corr_type='cf'):
        #Convert distances to alphas.
        at = DA_over_rd/self.DA_over_rd_fid
        ap = c_over_Hrd/self.c_over_Hrd_fid
        #With the new alphas, get the log likelihood.
        chi2 = self.interpolators[corr_type](at,ap)
        return chi2
