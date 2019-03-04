import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, RectBivariateSpline

basedir = '/Users/jfarr/Projects/H0LyA/'
scan_locations =   {'cf': basedir+'/data/BOSS_scans/BOSSDR11LyaF_k.scan',
                    'xcf': basedir+'/data/BOSS_scans/BOSSDR11QSOLyaF.scan'
                    }

BOSS_gaussian = {}
BOSS_gaussian_cf = {}
BOSS_gaussian_cf['z'] = 2.34
BOSS_gaussian_cf['par'] = 9.17
BOSS_gaussian_cf['par_sig'] = 0.28
BOSS_gaussian_cf['tra'] = 11.28
BOSS_gaussian_cf['tra_sig'] = 0.65
BOSS_gaussian_xcf = {}
BOSS_gaussian_xcf['z'] = 2.34
BOSS_gaussian_xcf['par'] = 9.0
BOSS_gaussian_xcf['par_sig'] = 0.3
BOSS_gaussian_xcf['tra'] = 10.8
BOSS_gaussian_xcf['tra_sig'] = 0.4
BOSS_gaussian['cf'] = BOSS_gaussian_cf
BOSS_gaussian['xcf'] = BOSS_gaussian_xcf

class chi2_interpolators():
    def __init__(self,scan_locations,DA_over_rd_fid,c_over_Hrd_fid,kind='linear'):
        #Create a dictionary containing an interpolator for each scan.
        interpolators = {}
        for corr_type in scan_locations:
            scan = np.loadtxt(scan_locations[corr_type])
            at = np.array(sorted(set(scan[:,0])))
            ap = np.array(sorted(set(scan[:,1])))
            grid = np.zeros((at.shape[0],ap.shape[0]))
            #for each ap value:
            for i in range(ap.shape[0]):
                #Filter the scan to only those corresponding to the ap value
                indices = scan[:,1]==ap[i]
                scan_chunk = scan[indices,:]
                #Ensure that they're sorted by at value
                scan_chunk = scan_chunk[scan_chunk[:,0].argsort()]
                #Add the chi2 column to the grid
                #Note the grid is of shape (at.shape,ap.shape)
                grid[:,i] = scan_chunk[:,2]

            interpolators[corr_type] = RectBivariateSpline(at,ap,grid,kx=1,ky=1)

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

<<<<<<< HEAD
    def plot_scan(self,corr_type,ap_min=0.8,ap_max=1.3,at_min=0.6,at_max=1.5,N_ap=100,N_at=100,sigma_contours=[1.,2.,3.],add_gaussian=True,savefig=False):

        #Get a chi grid of the desired properties
        at_grid = np.linspace(at_min,at_max,N_at)
        ap_grid = np.linspace(ap_min,ap_max,N_ap)
        chi2 = self.interpolators[corr_type](at_grid,ap_grid).T
        vmax = 5*(np.max(chi2//5) + 1)
        plt.imshow(chi2,aspect='auto',origin='lower',vmin=0.,vmax=vmax,extent=[at_min,at_max,ap_min,ap_max])
        plt.colorbar()

        #Add contours of the scan
        locs = [2.3,4.61,9.21]
        contour_locator = matplotlib.ticker.FixedLocator(locs)
        plt.contour(at_grid,ap_grid,chi2,locator=contour_locator,colors='white')

        #Add the same sigma contours of the gaussian approximation
        gaussian_ap = ((ap_grid - (BOSS_gaussian[corr_type]['par']/self.c_over_Hrd_fid))/(BOSS_gaussian[corr_type]['par_sig']/self.c_over_Hrd_fid))**2
        gaussian_at = ((at_grid - (BOSS_gaussian[corr_type]['tra']/self.DA_over_rd_fid))/(BOSS_gaussian[corr_type]['tra_sig']/self.DA_over_rd_fid))**2
        gaussian_grid = -np.log(np.outer(np.exp(-gaussian_ap),np.exp(-gaussian_at)))
        plt.contour(at_grid,ap_grid,gaussian_grid,locator=contour_locator,colors='white',linestyles='dotted')

        #Label axes, add a title, save and show
        plt.xlabel(r'$\alpha_t$')
        plt.ylabel(r'$\alpha_p$')
        plt.title('{} scan'.format(corr_type))
        if savefig:
            plt.savefig('{}_chi2_interpolation.pdf'.format(corr_type))
        plt.show()

        return
=======
class interpolate_2D():
    def __init__(self,data):
        self.data = data
        self.x = np.array(list(sorted(set(data[:,0]))))
        self.y = np.array(list(sorted(set(data[:,1]))))
        self.grid = np.zeros((self.y.shape,self.x.shape))
        #for each y value:
        for i in range(self.y.shape):
            #Filter the data to only those corresponding to the y value
            data_chunk = data[:,:][data[:,1]==self.y[i]]
            #Ensure that they're sorted by x value
            data_chunk = data_chunk[data_chunk[:,0].argsort()]
            #Add the chi2 column to the grid
            self.grid[i,:] = data_chunk[:,3]
        return
    def f_interp(self,x,y):
        i_upper = np.searchsorted(self.x,x)
        x1 = self.x[i_upper]
        x0 = self.x[i_upper - 1]
        j_upper = np.searchsorted(self.y,y)
        y1 = self.y[j_upper]
        y0 = self.y[j_upper - 1]
        y1_filter = self.data[:,1]==y1
        y0_filter = self.data[:,1]==y0
        fu = np.interp(x,self.data[:,0][y1_filter],self.data[:,2][y1_filter])
        fb = np.interp(x,self.data[:,0][y0_filter],self.data[:,2][y0_filter])
        f = fb + (fu - fb)*(y - y0)/(y1 - y0)
        return f

new=interpolate_2D(scan_cf)

for ap in ap_values:
    old_val = old.get_chi2_distances(10,ap*10)
    new_val = new.f_interp(1.0,ap)
    print(old_val,new_val)

import time
N=100000
start=time.time()
for i in range(N):
    old_val = interpolator(x_new,y_new)


print('interp2d took {:4.2f} seconds'.format(time.time()-start))
start=time.time()
for i in range(N):
    new_val = RBS(x_new,y_new)


print('RBS took {:4.2f} seconds'.format(time.time()-start))


def plot_scan_oldvsnew(old,new,N_at,N_ap):
    at_grid = np.linspace(0.6,1.5,N_at)
    ap_grid = np.linspace(0.8,1.3,N_ap)
    chi2_old = np.zeros((N_ap,N_at))
    chi2_new = np.zeros((N_ap,N_at))
    for i in range(N_ap):
        for j in range(N_at):
            chi2_old[i,j]=old.get_chi2_distances(at_grid[j]*da_over_rd_fid,ap_grid[i]*c_over_Hrd_fid)
            chi2_new[i,j]=new.f_interp(at_grid[j],ap_grid[i])
    plt.imshow(chi2_old,aspect='auto',origin='lower',vmin=0.,vmax=40.,extent=[min(at_grid),max(at_grid),min(ap_grid),max(ap_grid)])
    plt.colorbar()
    plt.savefig('old_chi2_interpolation_{}_{}.pdf'.format(N_at,N_ap))
    plt.show()
    plt.imshow(chi2_new,aspect='auto',origin='lower',vmin=0.,vmax=40.,extent=[min(at_grid),max(at_grid),min(ap_grid),max(ap_grid)])
    plt.colorbar()
    plt.savefig('new_chi2_interpolation_{}_{}.pdf'.format(N_at,N_ap))
    plt.show()
>>>>>>> 5f145757e9a2c955689ee3f75c83dc10387de86b
