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
