import numpy as np
import sys
import os
from cobaya.yaml import yaml_load_file
from cobaya import run

import random

model = '/Users/yg/chains/hede_mc/boxphi_f_mc_p20des.updated.yaml'
covmat = '/Users/yg/chains/hede_mc/boxphi_f_mc_p20des.covmat'
modelname = model.split(sep='/')[-1].split(sep='.')[0]

samples = ['/Users/yg/code/eftcamb_osc/boxphi_final/boxphi_f_mc_p20des_1.txt']

# minimizer = {'minimize':{'best_of':1, 'covmat': covmat}}
minimizer = {'minimize':{'best_of':1, 'covmat': covmat, 'override_scipy': None}}
# minimizer = {'minimize':{'best_of':1, 'covmat': covmat, 'override_iminuit': None}}
basepar = yaml_load_file(model)
basepar['sampler'] = minimizer
basepar['theory']['hede_shoot']['cosmo_path'] = '/Users/yg/code/eftcamb_osc'
basepar['theory']['hede_shoot']['python_path'] = '/Users/yg/code/eftcamb_osc'
basepar['theory']['camb']['path'] = '/Users/yg/code/eftcamb_osc'
basepar['likelihood']['pantheonplus']['data_path'] = '/Users/yg/data/PantheonPlus/'
basepar['likelihood']['pantheonplus']['python_path'] = '/Users/yg/code/eftcamb_osc'

minimum = 1e30
minimum_par = None

repeat = 3

idx = 1
for sp in samples:
    fin = open(sp, 'r')
    nms = fin.readline().replace('#','').split()
    fin.close()
    vec = np.loadtxt(sp)
    if 'xiv_phi' in nms:
        i=0
        while i < repeat:
            basepar['output'] = 'minimum/%s_%d'%( modelname, idx )
            for par in nms:
                if par in ['ombh2','omch2','H0']:
                    basepar['params'][par]['ref'] = vec[nms.index(par)]*random.normalvariate(1, 0.006)
                elif par in basepar['params'].keys():
                    basepar['params'][par]['ref'] = vec[nms.index(par)]
            updated_info, sampler = run(basepar, force=True)
            locmin = sampler.products()['minimum']
            if locmin['minuslogpost'] < minimum:
                minimum = locmin['minuslogpost']
                minimum_par = locmin.copy()
            i += 1
            idx += 1
    else:
        for val in [0.05, -0.05, -0.1]:
            basepar['output'] = 'minimum/%s_%d'%( modelname, idx )
            basepar['params']['xiv_phi']['ref'] = val
            for par in nms:
                if par in ['ombh2','omch2','H0']:
                    basepar['params'][par]['ref'] = vec[nms.index(par)]*random.normalvariate(1, 0.005)
                elif par in basepar['params'].keys():
                    basepar['params'][par]['ref'] = vec[nms.index(par)]
            updated_info, sampler = run(basepar, force=True)
            locmin = sampler.products()['minimum']
            if locmin['minuslogpost'] < minimum:
                minimum = locmin['minuslogpost']
                minimum_par = locmin.copy()
            idx += 1

print(minimum)
print(minimum_par)      