import numpy as np
from scipy.interpolate import RectBivariateSpline
from matplotlib import ticker
import matplotlib.pyplot as plt

#For reference.
basedir = '/Users/jfarr/Projects/H0LyA/'
scan_locations =   {'cf': basedir+'/data/BOSS_scans/BOSSDR11LyaF_k.scan',
                    'xcf': basedir+'/data/BOSS_scans/BOSSDR11QSOLyaF.scan'
                    }

#BOSS data Gaussian approximation. Could be put into a data file and stored
#separately.
BOSS_gauss = {}
BOSS_gauss_cf = {}
BOSS_gauss_cf['z'] = 2.34
BOSS_gauss_cf['par'] = 9.17
BOSS_gauss_cf['par_sig'] = 0.28
BOSS_gauss_cf['tra'] = 11.28
BOSS_gauss_cf['tra_sig'] = 0.65
BOSS_gauss_xcf = {}
BOSS_gauss_xcf['z'] = 2.34
BOSS_gauss_xcf['par'] = 9.0
BOSS_gauss_xcf['par_sig'] = 0.3
BOSS_gauss_xcf['tra'] = 10.8
BOSS_gauss_xcf['tra_sig'] = 0.4
BOSS_gauss['cf'] = BOSS_gauss_cf
BOSS_gauss['xcf'] = BOSS_gauss_xcf

#Class to read alpha_t by alpha_p chi2 scans e.g. from BOSS and interpolate.
class chi2_interpolators():
    def __init__(self,scan_locations,DA_over_rd_fid,c_over_Hrd_fid):
        """
        Arguments:
        scan_locations: dictionary of filepaths to the different scans, with
                        keys as scan types.
        DA_over_rd_fid: fiducial value of DA/rd used to calculate alpha_t.
        c_over_Hrd_fit: fiducial value of c/Hrd used to calculate alpha_p.
        """

        #Create a dictionary containing an interpolator for each scan.
        interpolators = {}
        for corr_type in scan_locations:
            scan = np.loadtxt(scan_locations[corr_type])

            at = np.array(sorted(set(scan[:,0])))
            ap = np.array(sorted(set(scan[:,1])))

            N_at = at.shape[0]
            N_ap = ap.shape[0]
            grid = np.zeros((N_at,N_ap))

            for i in range(N_ap):
                #Filter the data to only those corresponding to the ap value.
                indices = (scan[:,1]==ap[i])
                scan_chunk = scan[indices,:]
                #Ensure that they're sorted by at value.
                scan_chunk = scan_chunk[scan_chunk[:,0].argsort()]
                #Add the chi2 column to the grid.
                #Note that the grid is of shape (N_at,N_ap)
                grid[:,i] = scan_chunk[:,2]

            #Make the interpolator (x refers to at, y refers to ap).
            interpolators[corr_type] = RectBivariateSpline(at,ap,grid,kx=1,ky=1)

        #Add the dictionary to the object.
        self.interpolators = interpolators
        self.DA_over_rd_fid = DA_over_rd_fid
        self.c_over_Hrd_fid = c_over_Hrd_fid

        return

    #Function to return the interpolated value of chi2 given distance measures.
    def get_chi2_distances(self,DA_over_rd,c_over_Hrd,corr_type='cf'):
        """
        Arguments:
        DA_over_rd_fid: value of DA/rd to evaluate chi2 for.
        c_over_Hrd_fit: value of c/Hrd to evaluate chi2 for.
        corr_type:      which scan to interpolate.

        Returns:
        chi2:           value of chi2
        """

        #Convert distances to alphas.
        at = DA_over_rd/self.DA_over_rd_fid
        ap = c_over_Hrd/self.c_over_Hrd_fid

        #With the new alphas, get the log likelihood.
        chi2 = self.interpolators[corr_type](at,ap)

        return chi2

    #Function to plot the scan, interpolated to have the desired number of
    #pixels in each dimension. Can plot in terms of alpha or distances.
    def plot_interpolated_scan(self,N_at=100,N_ap=100,at_min=0.6,at_max=1.5,
                                ap_min=0.8,ap_max=1.3,mode='alphas',save=False,
                                corr_type='cf',add_sigma_contours=True,
                                add_gaussian_contours=True):
        """
        Arguments:
        N_at:                   Number of pixels in alpha_t (x direction).
        N_ap:                   Number of pixels in alpha_p (y direction).
        at_min:                 Min value of at to plot.
        at_max:                 Max value of at to plot.
        ap_min:                 Min value of ap to plot.
        ap_max:                 Max value of ap to plot.
        mode:                   Plot in terms of alpha ('alphas') or distance
                                ('distances')
        save:                   Would you like the plot to be saved?
        corr_type:              Which scan would you like to plot ('cf' or 'xcf'
                                normally).
        add_sigma_contours:     Would you like to add 1- 2- and 3- sigma
                                contours?
        add_gaussian_contours:  Would you like to add the same contours for the
                                Gaussian approximation?
        """

        #Genearte the at and ap values, then interpolate.
        at_grid = np.linspace(at_min,at_max,N_at)
        ap_grid = np.linspace(ap_min,ap_max,N_ap)
        interpolated_scan = self.get_chi2_distances(at_grid*self.DA_over_rd_fid,
                                                    ap_grid*self.c_over_Hrd_fid,
                                                    corr_type=corr_type).T

        #Make the figure.
        plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')

        #Determine what to put on the axes.
        if mode == 'alphas':
            x_grid = at_grid
            y_grid = ap_grid
            plt.xlabel(r'$\alpha_t$')
            plt.ylabel(r'$\alpha_p$')
        elif mode == 'distances':
            x_grid = at_grid * self.DA_over_rd_fid
            y_grid = ap_grid * self.c_over_Hrd_fid
            plt.xlabel(r'$D_A/r_d$')
            plt.ylabel(r'$c/Hr_d$')

        #Determine the top of the colourbar. Bottom is always 0 at the moment.
        max_chi2 = np.max(interpolated_scan)
        unit = 10. ** (np.ceil(np.log10(max_chi2))-1.)
        vmax = unit * np.ceil(max_chi2/unit)

        #Add 1- 2- and 3-sigma contours if desired.
        if add_sigma_contours:
            locs = [2.3,4.61,9.21]
            contour_locator = ticker.FixedLocator(locs)
            plt.contour(x_grid,y_grid,interpolated_scan,locator=contour_locator,
                        colors='white')

            #Add the same contours for the Gaussian log-likelihood if desired.
            if add_gaussian_contours:
                from scipy.stats import norm
                #Construct the Gaussians in each diection.
                ap_gauss_mean = BOSS_gauss[corr_type]['par']/self.c_over_Hrd_fid
                ap_gauss_sig = BOSS_gauss[corr_type]['par_sig']/self.c_over_Hrd_fid
                gauss_ap = norm.pdf(ap_grid,loc=ap_gauss_mean,scale=ap_gauss_sig)
                at_gauss_mean = BOSS_gauss[corr_type]['tra']/self.DA_over_rd_fid
                at_gauss_sig = BOSS_gauss[corr_type]['tra_sig']/self.DA_over_rd_fid
                gauss_at = norm.pdf(at_grid,loc=at_gauss_mean,scale=at_gauss_sig)

                #Combine and plot the contours.
                gauss_grid = -2.0 * np.log(2 * np.pi * ap_gauss_sig * at_gauss_sig *
                                        np.outer(gauss_ap,gauss_at))
                plt.contour(x_grid,y_grid,gauss_grid,locator=contour_locator,
                            colors='white',linestyles='dotted')

        #Show the plot, and save it if desired.
        plt.imshow(interpolated_scan,aspect='auto',origin='lower',vmin=0.,
                    vmax=vmax,extent=[min(x_grid),max(x_grid),min(y_grid),
                    max(y_grid)])
        plt.colorbar()
        if save:
            plt.savefig('interpolated_scan_{}_{}_{}.pdf'.format(N_at,N_ap,mode))
        plt.show()

        return

#Home-made interpolator class - basic but it works.
#RectBivariateSpline is used instead in the interpolators.
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
