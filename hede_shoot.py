from typing import Iterable, Sequence, Tuple, Union
from cobaya.typing import InfoDict
import numpy as np
from scipy.optimize import root

import pandas as pd

import os

from cobaya.theory import Theory
from cobaya.log import LoggedError
from cobaya.yaml import yaml_load_file
from cobaya.run import run
from cobaya.theories.cosmo.boltzmannbase import PowerSpectrumInterpolator
from cobaya.conventions import Const
from cobaya.likelihood import Likelihood

import sys

class hede_shoot(Theory):
    cosmo_path = None
    gravity_model = 3
    extra_args = {}

    def initialize(self):
        if self.cosmo_path != None:
            sys.path.insert(1, self.cosmo_path)
        import camb
        self.cosmo_shooter = camb

        model_param_num = 0
        use_ic = False
        evolve_hubble = False
        evolve_h = False
        bkint_pts = 10000
        if self.gravity_model == 3:
            print('==========================')
            print('Shooting initialized for model: XboxPhi EDE')
            model_param_num = 2
            use_ic = True
            bkint_pts = 3000
        else:
            raise Exception('Unkown model %s' % self.gravity_model)
        
        self.cosmo_pars = {'ombh2': None, 'omch2': None, 'H0': None}
        self.extargs = {'dark_energy_model': 'EFTCAMB',
                        'EFTflag' :4,
                        'FullMappingEFTmodel': 0,
                        'Horndeski_model': self.gravity_model, 
                        'Horndeski_parameter_number': model_param_num, 
                        'Horndeski_model_specific_ic': use_ic,
                        'Horndeski_evolve_hubble': evolve_hubble,
                        'EFT_evolve_delta_phi': True,
                        'EFT_evolve_metric_h': evolve_h,
                        'Horndeski_background_interpolation_num_points': bkint_pts,
                        'Horndeski_background_a_ini': 1e-14,
                        'EFTCAMB_turn_on_time': 1e-10,
                        'Horndeski_shooting': True,
                        'EFT_skip_stability': True,
                        }
        self.extargs.update(self.extra_args)
        self.f_hede = 0.1
        self.z_c = 5000
        self.model_par = -0.2

    def get_can_provide_params(self):
        if self.gravity_model == 3:
            ps = ["Horndeski_param1", "Horndeski_param2", "Honrdeski_phi_ini"]
        return ps
    
    def get_requirements(self):
        rqs = self.cosmo_pars.copy()
        rqs.update({'f_hede': None, 'lnzc': None})
        if self.gravity_model == 3:
            rqs['xiv_phi'] = None
        return rqs

    def calculate(self, state, want_derived=True, **params_values_dict):
        params = {}
        for par in self.cosmo_pars.keys():
            params[par] = params_values_dict[par]

        self.f_hede = params_values_dict['f_hede']
        self.z_c = np.exp(params_values_dict['lnzc'])
        
        if self.gravity_model == 3:
            self.model_par = params_values_dict['xiv_phi']

        x0 = np.empty(2)
        rho_fid = 1.21e-11*(self.z_c+1)**4 + 4.756e-8*(self.z_c+1)**3
        x0[0] = rho_fid/16/self.f_hede
        x0[1] = 2*np.sqrt(self.f_hede)
        # print("Asked for f=%.6e, zc=%.6e. Shooting started from v0=%.6e, phi_i=%.6e."%(self.f_hede, self.z_c, x0[0], x0[1]))

        try:
            rlt = root(self.is_ic, x0, tol=0.1, args=(params,), options={'eps': 0.05})
        except:
            return False
        if not rlt.success:
            self.log.debug("Shooting for horndeski parameters failed. "
                           "Assigning 0 likelihood and going on. ")
            return False
        v0 = rlt.x[0]
        phi_i = rlt.x[1]
        if self.gravity_model == 3:
            v_phi = 4.*v0*phi_i*phi_i*phi_i
            xi = self.model_par/v_phi

        # print("Shooting completed with v0=%.6e, phi_i=%.6e."%(v0, phi_i))

        if self.gravity_model == 3:
            state['derived'] = {'Horndeski_param1': xi, 'Horndeski_param2': v0, 'Honrdeski_phi_ini': phi_i}

        return True

    def is_ic(self, x, params):
        v0 = x[0]
        phi_i = x[1]
        if self.gravity_model == 3:
            v_phi = 4.*v0*phi_i*phi_i*phi_i
            xi = self.model_par/v_phi
        locparams = self.extargs.copy()
        locparams.update({'Horndeski_param1': xi, 'Horndeski_param2': v0, 'Honrdeski_phi_ini': phi_i})
        pars = self.cosmo_shooter.set_params(H0=params['H0'], ombh2=params['ombh2'], omch2=params['omch2'], **locparams)
        res = self.cosmo_shooter.get_background(pars, no_thermo=True)
        rlt = np.empty(2)
        rlt[0] = res.Params.EFTCAMB_parameter_cache.maxfrac_hdsk/self.f_hede - 1
        rlt[1] = res.Params.EFTCAMB_parameter_cache.z_maxfrac_hdsk/self.z_c - 1
        return rlt

