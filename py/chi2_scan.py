import numpy as np
from scipy.interpolate import interp2d, RectBivariateSpline

basedir = '/Users/jfarr/Projects/H0LyA/'
scan_locations =   {'cf': basedir+'/data/BOSS_scans/BOSSDR11LyaF_k.scan',
                    'xcf': basedir+'/data/BOSS_scans/BOSSDR11QSOLyaF.scan'
                    }

class chi2_interpolators():
    def __init__(self,scan_locations,DA_over_rd_fid,c_over_Hrd_fid,kind='linear'):
        #Create a dictionary containing an interpolator for each scan.
        interpolators = {}
        for corr_type in scan_locations:
            scan = np.loadtxt(scan_locations[corr_type])

            x = sorted(set(data[:,0]))
            y = sorted(set(data[:,1]))
            grid = np.zeros((y.shape,x.shape))
            #for each y value:
            for i in range(y.shape):
                #Filter the data to only those corresponding to the y value
                scan_chunk = scan[:,:][scan[:,1]==y[i]]
                #Ensure that they're sorted by x value
                scan_chunk = scan_chunk[scan_chunk[:,0].argsort()]
                #Add the chi2 column to the grid
                grid[i,:] = scan_chunk[:,3]
            interpolators[corr_type] = RectBivariateSpline(x,y,grid,kx=1,ky=1)

            #interpolators[corr_type] = interp2d(scan[:,0],scan[:,1],scan[:,2],kind='linear')
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
