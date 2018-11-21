import os
import numpy as np
from montepython.likelihood_class import Likelihood
import scipy.constants as const

class bao_lya_gaussian(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # define array for values of z and data points
        self.z = np.array([], 'float64')
        self.data = np.array([], 'float64')
        self.error = np.array([], 'float64')
        self.type = np.array([], 'int')

        # read redshifts and data points
        for line in open(os.path.join(
                self.data_directory, self.file), 'r'):
            if (line.strip().find('#') == -1) and (len(line.strip())>0):
                self.z = np.append(self.z, float(line.split()[0]))
                self.data = np.append(self.data, float(line.split()[1]))
                self.error = np.append(self.error, float(line.split()[2]))
                self.type = np.append(self.type, int(line.split()[3]))

        # number of data points
        self.num_points = np.shape(self.z)[0]

        # end of initialization

    # compute likelihood

    def loglkl(self, cosmo, data):

        chi2 = 0.

        # for each point, compute angular distance da, radial distance dr,
        # volume distance dv, sound horizon at baryon drag rs_d,
        # theoretical prediction and chi2 contribution
        # classes: (D_V/rs=3, Dv/Mpc=4, DA/rs=5, c/Hrs=6, rs/D_v=7, D_M/rs=8, H rs/rs_fid=9, D_M rs_fid/rs=10)
        for i in range(self.num_points):

            da = cosmo.angular_distance(self.z[i])
            dr = self.z[i] / cosmo.Hubble(self.z[i])
            #Don't understand the scaling here
            H  = cosmo.Hubble(self.z[i]) * const.c / 1000.
            #H  = cosmo.Hubble(self.z[i])

            dv = pow(da * da * (1 + self.z[i]) * (1 + self.z[i]) * dr, 1. / 3.)
            dm = da * (1 + self.z[i])

            rd = cosmo.rs_drag() * self.rd_rescale

            if (self.type[i] == 3):
                theo = dv / rd

            elif (self.type[i] == 4):
                theo = dv

            elif (self.type[i] == 5):
                theo = da / rd

            elif (self.type[i] == 6):
                theo = (const.c / 1000.) / (H * rd)

            elif (self.type[i] == 7):
                theo = rs / dv

            elif (self.type[i] == 8):
                theo = dm / rd

            elif (self.type[i] == 9):
                theo = (H * rd) / self.boss_rd_fid_in_Mpc

            elif (self.type[i] == 10):
                theo = (dm * self.boss_rd_fid_in_Mpc) / rd

            else:
                raise io_mp.LikelihoodError(
                    "In likelihood %s. " % self.name +
                    "BAO data type %s " % self.type[i] +
                    "in %d-th line not understood" % i)

            chi2 += ((theo - self.data[i]) / self.error[i]) ** 2

        # return ln(L)
        lkl = - 0.5 * chi2

        return lkl
