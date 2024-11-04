import numpy as np

import pandas as pd

import os

from cobaya.likelihood import Likelihood

class pantheonplus(Likelihood):
    data_path = "/data2/data/PantheonPlus/"
    with_sh0es = False

    def initialize(self):
        print("== Setting up PantheonPlus Likelihood ==")

        filename = os.path.join(self.data_path, 'Pantheon+SH0ES.dat')
        print("Loading data from {}".format(filename))
        data = pd.read_csv(filename, sep='\s+')
        self.origlen = len(data)

        if self.with_sh0es:
            self.ww = (data['zHD']>0.01) | (np.array(data['IS_CALIBRATOR'],dtype=bool))
            self.is_calibrator = np.array(data['IS_CALIBRATOR'][self.ww], dtype=bool)
            self.cepheid_distance = data['CEPH_DIST'][self.ww]
        else:
            self.ww = (data['zHD'] > 0.01)

        # use the vpec corrected redshift for zCMB
        self.zCMB = data['zHD'][self.ww]
        self.zHEL = data['zHEL'][self.ww]
        self.m_obs = data['m_b_corr'][self.ww]

        self.invcov = self._build_invcov()

        print('== PantheonPlus likelihood setup done ==')

    def get_requirements(self):
        # State requisites to the theory code
        reqs = {"angular_diameter_distance": {"z": self.zCMB}}
        reqs["Mb"] = None
        return reqs

    def logp(self, **params_values):
        da_th = self.provider.get_angular_diameter_distance(self.zCMB)
        m_th = 5.0*np.log10((1.0+self.zCMB)*(1.0+self.zHEL)*da_th)+25.
        if self.with_sh0es:
            m_th[self.is_calibrator] = self.cepheid_distance[self.is_calibrator]
        Mb = params_values.get('Mb', None)
        m_th += Mb
        diffmag = m_th - self.m_obs
        return -0.5*np.dot(diffmag, np.dot(self.invcov, diffmag))

    def _build_invcov(self):
        filename = os.path.join(self.data_path, 'Pantheon+SH0ES_STAT+SYS.cov')
        print("Loading covariance from {}".format(filename))

        # The file format for the covariance has the first line as an integer
        # indicating the number of covariance elements, and the the subsequent
        # lines being the elements.
        # This function reads in the file and the nasty for loops trim down the covariance
        # to match the only rows of data that are used for cosmology

        f = open(filename)
        line = f.readline()
        n = int(len(self.zCMB))
        cov = np.zeros((n, n))
        ii = -1
        jj = -1
        mine = 999
        maxe = -999
        for i in range(self.origlen):
            jj = -1
            if self.ww[i]:
                ii += 1
            for j in range(self.origlen):
                if self.ww[j]:
                    jj += 1
                val = float(f.readline())
                if self.ww[i]:
                    if self.ww[j]:
                        cov[ii, jj] = val
        f.close()

        return np.linalg.inv(cov)
