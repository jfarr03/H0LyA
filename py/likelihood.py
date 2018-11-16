import numpy as np
from scipy.interpolate import interp2d

basedir = '../BOSS_data/'
data_locations =   {'cf': basedir+'/BOSSDR11LyaF_k.scan'
                    'xcf': basedir+'/BOSSDR11QSOLyaF.scan'
                    }

def get_loglikelihood_distances():

    return
    
def get_loglikelihood_alphas(ap,at,corr_type):

    #Load the data.
    data = np.loadtxt(data_locations[corr_type])

    # TODO: which way round are ap and at? Does it matter
    data_ap = data[:,0]
    data_at = data[:,1]
    data_chi2 = data[:,2]

    interpolator = interp2d(data_ap,data_at,data_chi2)
    chi2 = interpolator(ap,at)

    return -chi2/2.

def distances_to_alphas(DA_over_rd,c_over_rdH):

    #Calculate the fiducial values of the distances.

    #Convert to ratios (alphas).
    ap = DA_over_rd/DA_over_rd_fid
    at = c_over_rdH/c_over_rdH_fid

    return ap, at

def
