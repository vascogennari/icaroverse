import os, sys, configparser, shutil, time, re
from optparse import OptionParser
from inspect import getmembers, isclass
import multiprocessing as mp

import pickle, h5py, pandas as pd, json
import icarogw, bilby, astropy

if icarogw.cupy_pal.is_there_cupy(): import cupy  as xp
else:                                import numpy as xp

from numpy import minimum

# Internal imports
from . import options
from . import postprocessing



def get_wrapper(wrap_name, input_wrapper = None, order = None, transition = None, smoothing = None, z_mixture = None, cosmo_wrap = False, bkg_cosmo_wrap_name = None, n_splines = None, spacing = None, zmax = 20.):

    print('\t{}'.format(wrap_name))
    wrap = getattr(icarogw.wrappers, wrap_name)
    if cosmo_wrap:
        # if bkg_cosmo_wrap_name is not None, it is assumed that wrap_name refers to a modified gravity cosmology wrapper
        if bkg_cosmo_wrap_name is not None:
            bkg_wrap = get_wrapper(bkg_cosmo_wrap_name, cosmo_wrap=True, zmax=zmax)
            return wrap(bkg_wrap)
        else:
            return wrap(zmax=zmax)
    elif transition == None:
        if order == None:
            if not input_wrapper == None:

                return wrap(input_wrapper)
            elif wrap_name == 'PowerLaw_PowerLaw' or wrap_name == 'PowerLaw_PowerLaw_PowerLaw' or wrap_name == 'PowerLaw_PowerLaw_PowerLaw_PowerLaw' or wrap_name == 'PowerLaw_PowerLaw_Gaussian':
                return wrap(flag_powerlaw_smoothing = smoothing)
            elif 'Spline' in wrap_name:
                print('\t\tUsing a spline model with {} basis elements. Knots spacing: {}.\n'.format(n_splines, spacing))
                return wrap(n_basis = n_splines, spacing = spacing)
            else:
                return wrap()
        else:
            # GaussianRedshift-order-X model.
            return wrap(order = order)
    else:
        if   wrap_name == 'PowerLaw_GaussianRedshiftLinear' or wrap_name == 'PowerLaw_GaussianRedshiftQuadratic' or wrap_name == 'PowerLaw_GaussianRedshiftPowerLaw' or wrap_name == 'PowerLaw_GaussianRedshiftSigmoid' or wrap_name == 'PowerLawBroken_GaussianRedshiftLinear' or wrap_name == 'PowerLawRedshiftLinear_GaussianRedshiftLinear' or wrap_name == 'PowerLaw_GaussianRedshiftLinear_GaussianRedshiftLinear' or wrap_name == 'PowerLawRedshiftLinear_PowerLawRedshiftLinear_PowerLawRedshiftLinear_PowerLawRedshiftLinear':
            return wrap(redshift_transition = transition, flag_redshift_mixture = z_mixture, flag_powerlaw_smoothing = smoothing)
        elif wrap_name == 'GaussianRedshiftLinear_GaussianRedshiftLinear' or wrap_name == 'GaussianRedshiftLinear_GaussianRedshiftLinear_GaussianRedshiftLinear' or wrap_name == 'PowerLawRedshiftLinear_PowerLawRedshiftLinear_PowerLawRedshiftLinear' or wrap_name == 'PowerLawRedshiftLinear_PowerLawRedshiftLinear_GaussianRedshiftLinear':
            return wrap(redshift_transition = transition, flag_redshift_mixture = z_mixture)
        elif wrap_name == 'DoublePowerlawRedshift':
            return wrap(redshift_transition = transition)

def print_dictionary(dictionary):
      
    for key in dictionary.keys():
        max_len = len(max(dictionary.keys(), key = len))
        if not key == 'all-priors': print('\t{}  {}'.format(key.ljust(max_len), dictionary[key]))

def check_effective_number_injections(pars, likelihood, n_events, maxL_values = None):

    def single_likelihood_eval():
        count = time.time()
        _     = likelihood.log_likelihood()
        print('\n\tA single likelihood evaluation takes {0:.5f} [s].'.format(time.time() - count))

    reference_model_dict = None
    if maxL_values == None:
        # if (not pars['real-data']) and (not pars['true-values'] == {}):
        if (not pars['true-values'] == {}):
            reference_model_dict = pars['true-values']
            reference_model_str  = 'injected'
    else:
        reference_model_dict = maxL_values
        reference_model_str  = 'maximum likelihood'

    if not reference_model_dict == None:
        # Set rate model parameters at true values 
        likelihood.parameters = {key: reference_model_dict[key] for key in likelihood.rate_model.population_parameters}
        # First likelihood evaluation at true values
        single_likelihood_eval()

        if pars['loglike-var'] is None or pars['loglike-var'] <= 0.:
            # Check effective number of injections.
            N_eff_inj = likelihood.injections.effective_injections_number()
            stability = N_eff_inj / (4 * n_events)
            print('\n\tThe effective number of injections for the {2} model is {0:.1f}. N_eff_inj/4*N_events is {1:.1f}.'.format(N_eff_inj, stability, reference_model_str))
            if stability < 1: print('\n\tWARNING: The number of injections is not enough to ensure numerical stability in the computation of selection effects in the likelihood. Please consider using a larger set of injections.')
            # Check effective numer of posterior samples.
            try:
                N_eff_PE  = xp.min(likelihood.posterior_samples_dict.get_effective_number_of_PE())
                print('\n\tThe minimum effective number of PE samples for the {1} model is {0:.1f}.'.format(N_eff_PE, reference_model_str))
            except AttributeError as err:
                # The first likelihood evaluation at true population values gives 0 because the effective number of injections is below threshold.
                # Consequently the initialisation of posterior samples weights is skipped.
                raise AttributeError(err, "* The effective number of injections for the true population values is below threshold.")

        else:
            try:
                N_eff_inj = likelihood.injections.effective_injections_number()
                stability = N_eff_inj / (4 * n_events)
                print('\n\tThe effective number of injections for the {2} model is {0:.1f}. N_eff_inj/4*N_events is {1:.1f}.'.format(N_eff_inj, stability, reference_model_str))
            except:
                print('\n\tWARNING: Could not retrieve Neffinj')
            try:
                N_eff_PEs = likelihood.posterior_samples_dict.get_effective_number_of_PE()
                N_eff_PE = xp.min(N_eff_PEs)
                print('\n\tThe minimum effective number of PE samples for the {1} model is {0:.1f}.'.format(N_eff_PE, reference_model_str))
            except:
                print('\n\tWARNING: Could not retrieve NeffPE')

            loglike_var = likelihood.likelihood_variance
            print('\n\tThe variance on the log-likelihood for the {1} model is {0:.2e}.'.format(loglike_var, reference_model_str))
            if likelihood.likelihood_variance_thr is not None and likelihood.likelihood_variance > likelihood.likelihood_variance_thr:
                raise ValueError("* The variance of the log likelihood for the true population values is below threshold.")


class Wrappers:

    def __init__(self, pars):

        self.pars = pars
        print('\n * Initialising the wrappers.\n')

    def PrimaryMass(self, pars):

        mp, ms = pars['model-primary'], pars['model-secondary']
        single_mass, smoothing, z_transition, z_mixture = pars['single-mass'], pars['low-smoothing'], pars['redshift-transition'], pars['redshift-mixture']
        # This is subject to be completed in the future with the addition of other primary mass distributions models to icarogw
        models = {
            'PowerLaw':                                                                                    {'wrap name': 'massprior_PowerLaw',                                                                          'z evolution': False, 'smoothing': 'global'},
            'PowerLaw-Gaussian':                                                                           {'wrap name': 'massprior_PowerLawPeak',                                                                      'z evolution': False, 'smoothing': 'global'},
            'PowerLaw-Gaussian-Gaussian':                                                                  {'wrap name': 'massprior_MultiPeak',                                                                         'z evolution': False, 'smoothing': 'global'},
            'PowerLaw-PowerLaw':                                                                           {'wrap name': 'PowerLaw_PowerLaw',                                                                           'z evolution': False, 'smoothing': 'component-wise'},
            'PowerLaw-PowerLaw-PowerLaw':                                                                  {'wrap name': 'PowerLaw_PowerLaw_PowerLaw',                                                                  'z evolution': False, 'smoothing': 'component-wise'},
            'PowerLaw-PowerLaw-PowerLaw-PowerLaw':                                                         {'wrap name': 'PowerLaw_PowerLaw_PowerLaw_PowerLaw',                                                         'z evolution': False, 'smoothing': 'component-wise'},
            'PowerLaw-PowerLaw-Gaussian':                                                                  {'wrap name': 'PowerLaw_PowerLaw_Gaussian',                                                                  'z evolution': False, 'smoothing': 'component-wise'},
            'PowerLaw-GaussianRedshiftLinear':                                                             {'wrap name': 'PowerLaw_GaussianRedshiftLinear',                                                             'z evolution': True,  'smoothing': 'component-wise'},
            'PowerLaw-GaussianRedshiftQuadratic':                                                          {'wrap name': 'PowerLaw_GaussianRedshiftQuadratic',                                                          'z evolution': True,  'smoothing': 'component-wise'},
            'PowerLaw-GaussianRedshiftPowerLaw':                                                           {'wrap name': 'PowerLaw_GaussianRedshiftPowerLaw',                                                           'z evolution': True,  'smoothing': 'component-wise'},
            'PowerLaw-GaussianRedshiftSigmoid':                                                            {'wrap name': 'PowerLaw_GaussianRedshiftSigmoid',                                                            'z evolution': True,  'smoothing': 'component-wise'},
            'PowerLawBroken-GaussianRedshiftLinear':                                                       {'wrap name': 'PowerLawBroken_GaussianRedshiftLinear',                                                       'z evolution': True,  'smoothing': 'component-wise'},
            'PowerLawRedshiftLinear-GaussianRedshiftLinear':                                               {'wrap name': 'PowerLawRedshiftLinear_GaussianRedshiftLinear',                                               'z evolution': True,  'smoothing': 'component-wise'},
            'PowerLaw-GaussianRedshiftLinear-GaussianRedshiftLinear':                                      {'wrap name': 'PowerLaw_GaussianRedshiftLinear_GaussianRedshiftLinear',                                      'z evolution': True,  'smoothing': 'component-wise'},
            'GaussianRedshiftLinear-GaussianRedshiftLinear':                                               {'wrap name': 'GaussianRedshiftLinear_GaussianRedshiftLinear',                                               'z evolution': True,  'smoothing': 'component-wise'},
            'GaussianRedshiftLinear-GaussianRedshiftLinear-GaussianRedshiftLinear':                        {'wrap name': 'GaussianRedshiftLinear_GaussianRedshiftLinear_GaussianRedshiftLinear',                        'z evolution': True,  'smoothing': 'component-wise'},
            'PowerLawRedshiftLinear-PowerLawRedshiftLinear-PowerLawRedshiftLinear':                        {'wrap name': 'PowerLawRedshiftLinear_PowerLawRedshiftLinear_PowerLawRedshiftLinear',                        'z evolution': True,  'smoothing': 'component-wise'},
            'PowerLawRedshiftLinear-PowerLawRedshiftLinear-PowerLawRedshiftLinear-PowerLawRedshiftLinear': {'wrap name': 'PowerLawRedshiftLinear_PowerLawRedshiftLinear_PowerLawRedshiftLinear_PowerLawRedshiftLinear', 'z evolution': True,  'smoothing': 'component-wise'},
            'PowerLawRedshiftLinear-PowerLawRedshiftLinear-GaussianRedshiftLinear':                        {'wrap name': 'PowerLawRedshiftLinear_PowerLawRedshiftLinear_GaussianRedshiftLinear',                        'z evolution': True,  'smoothing': 'component-wise'},
            'Gaussian':                                                                                    {'wrap name': 'Gaussian',                                                                                    'z evolution': False, 'smoothing': 'included'},
            'GaussianRedshift-order-X':                                                                    {'wrap name': 'GaussianEvolving',                                                                            'z evolution': True},
            'Splines-Quadratic':                                                                           {'wrap name': 'QuadraticSpline',                                                                             'z evolution': False, 'smoothing': 'component-wise'},
            'Splines-Cubic':                                                                               {'wrap name': 'CubicSpline',                                                                                 'z evolution': False, 'smoothing': 'component-wise'},
            'Uniform':                                                                                     {'wrap name': 'Uniform',                                                                                     'z evolution': False, 'smoothing': 'included'},
            'DoublePowerlaw':                                                                              {'wrap name': 'DoublePowerlaw',                                                                              'z evolution': False, 'smoothing': 'included'},
            'DoublePowerlaw-Gaussian':                                                                     {'wrap name': 'DoublePowerlaw_Gaussian',                                                                     'z evolution': False, 'smoothing': 'included'},
            'DoublePowerlawRedshift':                                                                      {'wrap name': 'DoublePowerlawRedshift',                                                                      'z evolution': True,  'smoothing': 'included'},
            'Johnson':                                                                                     {'wrap name': 'Johnson',                                                                                     'z evolution': False, 'smoothing': 'included'},            
        }
        # This is to make sure one can only use the models that are present in one's currently installed version of icarogw, AND that the present pipeline can handle.
        available_icarogw_models = dict(getmembers(icarogw.wrappers, isclass))
        icarogw_models = [m for m in models if models[m]['wrap name'] in available_icarogw_models]

        order = -1
        search = re.search("GaussianRedshift-order-(?P<order>[0-9]+)", mp)
        if search is not None: mp, order = 'GaussianRedshift-order-X', int(search.group('order'))

        if mp in icarogw_models:
            # Non-evolving models.
            if   (not models[mp]['z evolution']) and models[mp]['smoothing'] == 'global':
                w = get_wrapper(models[mp]['wrap name'])
                if (not (single_mass and 'Mass2' in ms)) and smoothing: w = get_wrapper('lowSmoothedwrapper', input_wrapper = w)
            elif 'Splines' in mp:
                w = get_wrapper(models[mp]['wrap name'], n_splines = pars['splines-number'], spacing = pars['spacing'])
            elif (not models[mp]['z evolution']) and models[mp]['smoothing'] == 'component-wise':
                w = get_wrapper(models[mp]['wrap name'], smoothing = smoothing)
            elif (not models[mp]['z evolution']) and models[mp]['smoothing'] == 'included':
                w = get_wrapper(models[mp]['wrap name'])
            # Evolving models.
            elif (    models[mp]['z evolution']) and order > 0: # GaussianRedshift-order-X model.
                w = get_wrapper(models[mp]['wrap name'],                        order = order,                                                 )
            elif (    models[mp]['z evolution']) and models[mp]['smoothing'] == 'included':
                w = get_wrapper(models[mp]['wrap name'],                                       transition = z_transition                       )
            elif (    models[mp]['z evolution']):
                w = get_wrapper(models[mp]['wrap name'], smoothing = smoothing,                transition = z_transition, z_mixture = z_mixture)
        # Unknown model
        else:
            raise ValueError("Unknown model for the Primary Mass distribution: {}.\nPlease choose from the available models:\n\t{}".format(mp, "\n\t".join(icarogw_models)))

        return w

    def SecondaryMass(self, pars, m1w = None):

        ms = pars['model-secondary']
        single_mass, smoothing = pars['single-mass'], pars['low-smoothing']
        # This is subject to be completed in the future with the addition of other primary mass distributions models to icarogw
        models = {
            'Mass2-PowerLaw':     {'wrap name': 'm1m2_conditioned',          'var': 'm2'},
            'MassRatio-Gaussian': {'wrap name': 'mass_ratio_prior_Gaussian', 'var': 'q' }, 
            'MassRatio-PowerLaw': {'wrap name': 'mass_ratio_prior_Powerlaw', 'var': 'q' }, 
            'MassRatio-Gamma':    {'wrap name': 'Gamma',                     'var': 'q' }, 
            'MassRatio-Beta':     {'wrap name': 'Beta',                      'var': 'q' },
        }
        # This is to make sure one can only use the models that are present in one's currently installed version of icarogw, AND that the present pipeline can handle.
        available_icarogw_models = dict(getmembers(icarogw.wrappers, isclass))
        icarogw_models = [m for m in models if models[m]['wrap name'] in available_icarogw_models]

        if single_mass:
            print('\t * Skipping secondary mass wrapper')
            w = None
        elif ms in icarogw_models:
            if   models[ms]['var'] == 'm2': w = get_wrapper(models[ms]['wrap name'] + '_lowpass'*smoothing, input_wrapper = m1w)
            elif models[ms]['var'] == 'q':  w = get_wrapper(models[ms]['wrap name']                                            )
        else:
            raise ValueError("Unknown model for the Secondary Mass distribution: {}.\nPlease choose from the available models:\n\t{}".format(ms, "\n\t".join(icarogw_models)))
        return w

    def RateEvolution(self, pars):

        mr = pars['model-rate']
        # This is subject to be completed in the future with the addition of other rate evolution models to icarogw
        models = {
            'MadauDickinson':                   {'wrap name': 'rateevolution_Madau'},
            'BetaDistribution':                 {'wrap name': 'rateevolution_beta'},
            'BetaDistribution-Line':            {'wrap name': 'rateevolution_beta_line'},
            'MadauDickinson-GammaDistribution': {'wrap name': 'rateevolution_Madau_gamma'},
            'PowerLaw':                         {'wrap name': 'rateevolution_PowerLaw'},
            'Gaussian':                         {'wrap name': 'rateevolution_Gaussian'},
            'RedshiftProbability-Beta':         {'wrap name': 'rateevolution_beta_redshift_probability'},
            'RedshiftProbability-Uniform':      {'wrap name': 'rateevolution_uniform_redshift_probability'},
            'RedshiftProbability-PowerLaw':     {'wrap name': 'rateevolution_powerlaw_redshift_probability'},
            'LuminosityProbability-Beta':       {'wrap name': 'rateevolution_beta_redshift_probability'},
        }
        # This is to make sure one can only use the models that are present in one's currently installed version of icarogw, AND that the present pipeline can handle.
        available_icarogw_models = dict(getmembers(icarogw.wrappers, isclass))
        icarogw_models = [m for m in models if models[m]['wrap name'] in available_icarogw_models]

        if mr in icarogw_models:
            w = get_wrapper(models[mr]['wrap name'])
        else:
            raise ValueError("Unknown model for the Rate Evolution: {}.\nPlease choose from the available models:\n\t{}".format(mr, "\n\t".join(icarogw_models)))
        return w
    
    def Cosmology(self, pars):

        mc, mb = pars['model-cosmology'], pars['model-bkg-cosmo']
        # This is subject to be completed in the future with the addition of other cosmological models to icarogw
        models = {
            'FlatLambdaCDM': {'wrap name': 'FlatLambdaCDM_wrap', 'class': 'GR', }, 
            'FlatwCDM':      {'wrap name': 'FlatwCDM_wrap',      'class': 'GR', }, 
            'Flatw0waCDM':   {'wrap name': 'Flatw0waCDM_wrap',   'class': 'GR', }, 
            'wIDS_linDE':    {'wrap name': 'wIDS_linDE_wrap',    'class': 'GR', }, 
            'Xi0':           {'wrap name': 'Xi0_mod_wrap',       'class': 'MG', }, 
            'eps0':          {'wrap name': 'eps0_mod_wrap',      'class': 'MG', }, 
            'extraD':        {'wrap name': 'extraD_mod_wrap',    'class': 'MG', }, 
            'cM':            {'wrap name': 'cM_mod_wrap',        'class': 'MG', }, 
            'alphalog':      {'wrap name': 'alphalog_mod_wrap',  'class': 'MG', }, 
        }
        # This is to make sure one can only use the models that are present in one's currently installed version of icarogw, AND that the present pipeline can handle.
        available_icarogw_models = dict(getmembers(icarogw.wrappers, isclass))
        icarogw_models = [m for m in models if models[m]['wrap name'] in available_icarogw_models]

        if   (mc in icarogw_models and models[mc]['class'] == 'MG') and (mb in icarogw_models and models[mc]['class'] == 'GR'):
            raise ValueError("Unknown GR model for the background cosmology: {}.\nPlease choose from the available models:\n\t{}".format(mb, "\n\t".join([m for m in icarogw_models if models[m]['class'] == 'GR'])))
        elif (mc in icarogw_models and models[mc]['class'] == 'MG'):
            w = get_wrapper(models[mc]['wrap name'], cosmo_wrap=True, bkg_cosmo_wrap_name=models[mb]['wrap name'], zmax=pars['zmax'])
        elif (mc in icarogw_models and models[mc]['class'] == 'GR'):
            w = get_wrapper(models[mc]['wrap name'], cosmo_wrap=True, bkg_cosmo_wrap_name=None                   , zmax=pars['zmax'])
        else:
            raise ValueError("Unknown model for the Cosmology: {}.\nPlease choose from the available models:\n\t{}".format(mc, "\n\t".join(icarogw_models)))
        return w

    def ReferenceCosmology(self):
        '''
            Reference cosmology used to generate the injections to compute selection effects.
            This wrapper is only used to convert the injections' prior and in the post-processing.
        '''
        w = icarogw.cosmology.astropycosmology(self.pars['zmax'])
        w.build_cosmology(astropy.cosmology.FlatLambdaCDM(H0 = self.pars['ref-cosmology']['H0'], Om0 = self.pars['ref-cosmology']['Om0']))
        return w

    def return_Wrappers(self):

        self.Wrapper_PrimaryMass   = self.PrimaryMass(self.pars)
        self.Wrapper_SecondaryMass = self.SecondaryMass(self.pars, self.Wrapper_PrimaryMass)
        self.Wrapper_RateEvolution = self.RateEvolution(self.pars)
        self.Wrapper_Cosmology     = self.Cosmology(self.pars)
        self.Wrapper_RefCosmology  = self.ReferenceCosmology()

        return self.Wrapper_PrimaryMass, self.Wrapper_SecondaryMass, self.Wrapper_RateEvolution, self.Wrapper_Cosmology, self.Wrapper_RefCosmology



class Rate():
      
    def __init__(self, pars, m1w, m2w, rw, cw):

        if not pars['single-mass']:
            if   not 'Redshift' in pars['model-primary'] and not 'MassRatio'  in pars['model-secondary']:
                self.w = icarogw.rates.CBC_vanilla_rate(                  cw,      m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('CBC_vanilla_rate'))
            elif not 'Redshift' in pars['model-primary'] and      'Gamma'     in pars['model-secondary'] and not 'Probability' in pars['model-rate']:
                self.w = icarogw.rates.MBH_rate(                          cw, m1w, m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('MBH_rate'))
            elif not 'Redshift' in pars['model-primary'] and 'Probability' in pars['model-rate']:
                self.w = icarogw.rates.MBH_redshift_rate(                 cw, m1w, m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('MBH_redshift_rate'))
            elif 'Probability' in pars['model-rate']:
                self.w = icarogw.rates.MBH_redshift_rate_given_redshift(  cw, m1w, m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('MBH_redshift_rate_given_redshift'))
            elif not 'Redshift' in pars['model-primary'] and      'MassRatio' in pars['model-secondary']:
                self.w = icarogw.rates.CBC_rate_m1_q(                     cw, m1w, m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('CBC_rate_m1_q'))
            elif     'Redshift' in pars['model-primary'] and  not 'MassRatio' in pars['model-secondary']:
                self.w = icarogw.rates.CBC_rate_m1_given_redshift_m2(     cw, m1w, m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('CBC_rate_m1_given_redshift_m2'))
            elif     'Redshift' in pars['model-primary'] and      'MassRatio' in pars['model-secondary']:
                self.w = icarogw.rates.CBC_rate_m1_given_redshift_q(      cw, m1w, m2w, rw, scale_free = pars['scale-free'])
                print('\t{}'.format('CBC_rate_m1_given_redshift_q'))
        else:
            if not 'Probability' in pars['model-rate']:
                self.w = icarogw.rates.CBC_rate_m_given_redshift(         cw, m1w,      rw, scale_free = pars['scale-free'])
                print('\t{}'.format('CBC_rate_m_given_redshift'))
            else:
                if not 'Luminosity' in pars['model-rate']:
                    self.w = icarogw.rates.CBC_redshift_rate_m_given_redshift(cw, m1w, rw, scale_free = pars['scale-free'])
                    print('\t{}'.format('CBC_redshift_rate_m_given_redshift'))
                else:
                    self.w = icarogw.rates.CBC_redshift_rate_m_given_luminosity(cw, m1w, rw)
                    print('\t{}'.format('CBC_redshift_rate_m_given_luminosity'))

        print('\n * Population parameters.\n')
        print('\t{}'.format('[%s]' % ', '.join(map(str, self.w.population_parameters))))
        pars['population-parameters'] = self.w.population_parameters

    def return_Rate(self):
        return self.w



class SelectionEffects:
        
    def __init__(self, pars, ref_cosmo):

        print('\n * Loading injections for selection effects.')

        # Use IGWN real-noise injections from IGWN.
        if pars['real-noise-injections']:

            print('\n\tUsing IGWN sensitivity estimates in real noise to evaluate selection effects.')

            from IGWN_pointers import sensitivity_estimates
            path = sensitivity_estimates[pars['catalog']].replace("~IGWN_injections_path", pars['injections-path'])
            try: data_inj = h5py.File(path)
            except: raise ValueError('Could not open the file containing the injections for selection effects. Please verify that you have downloaded the IGWN sensitivity estimates and that the path is correct:\n{}'.format(sensitivity_estimates[pars['catalog']]))
            print('\n\t{}'.format(path))

            # FIXME: Add control to make sure that the injections used correspond to the correct catalog.
            if   pars['catalog'] == 'GWTC-3':   obs_time = None # FIXME: Where can this value be found?
            elif pars['catalog'] == 'GWTC-4.0': obs_time = data_inj.attrs['total_analysis_time'] / 31557600.0
            elif pars['catalog'] == 'O3':       obs_time = (28519200 / 86400) / 365 # FIXME: Where is this read from and why is it hardcoded?
            elif pars['catalog'] == 'O4a':      obs_time = None # FIXME: Where can this value be found?
            else:
                raise ValueError('Unknown catalog option. Please choose from GWTC-3, GWTC-4.0, O3 or O4a.')

            pars['injections-number'] = data_inj.attrs['total_generated']

            if   pars['catalog'] == 'O3':

                prior  = icarogw.cupy_pal.np2cp(data_inj['injections/mass1_source_mass2_source_sampling_pdf'][()] * data_inj['injections/redshift_sampling_pdf'][()])
                # Converting the injections from source to detector frame, we need to correct the injections prior by the Jacobian of the transformation (m1s,m2s,z)->(m1d,m2d,dL).
                prior *= icarogw.conversions.source2detector_jacobian(icarogw.cupy_pal.np2cp(data_inj['injections/redshift'][()]), ref_cosmo)

                tmp = xp.vstack([data_inj['injections'][key] for key in ['ifar_cwb', 'ifar_gstlal', 'ifar_mbta', 'ifar_pycbc_bbh', 'ifar_pycbc_hyperbank']])
                ifarmax = xp.max(tmp, axis = 0)

                inj_dict = {
                    'mass_1':              xp.array(data_inj['injections/mass1'][()]),
                    'mass_2':              xp.array(data_inj['injections/mass2'][()]),
                    'luminosity_distance': xp.array(data_inj['injections/distance'][()])}

            elif pars['catalog'] == 'GWTC-4.0':

                pars['injections-number'] = data_inj.attrs['total_generated']
                events = data_inj['events'][:]
                ifarmax = 1 / xp.min(xp.array([events[search + '_far'] for search in data_inj.attrs['searches']]), axis = 0)
                if   pars['selection-effects-cut'] == 'snr' : real_noise_condition = ifarmax >= pars['snr-cut']
                elif pars['selection-effects-cut'] == 'ifar': real_noise_condition = ifarmax >= pars['ifar-cut']
                else:
                    raise ValueError('Unknown option to compute the selection effects cut.')
                selected_filt = (real_noise_condition) | (xp.array(events['semianalytic_observed_phase_maximized_snr_net']) >= pars['snr-cut-analytic'])

                lnprior = xp.array(events['lnpdraw_mass1_source_mass2_source_redshift_spin1_magnitude_spin1_polar_angle_spin1_azimuthal_angle_spin2_magnitude_spin2_polar_angle_spin2_azimuthal_angle'])
                lnprior -= xp.log(xp.sin(xp.array(events['spin2_polar_angle'])) * xp.sin(xp.array(events['spin1_polar_angle']))) # Accounts from implied Jacobian from t->cost.
                lnprior += xp.log(4 * xp.pi * xp.pi) # Azimuthal angles, not accouted in icarogw have uniform prior 1/2pi that we remove.
                lnprior -= xp.log(xp.array(events['weights'])) # Weights for different observing runs.
                lnprior -= xp.log(xp.array(events['dluminosity_distance_dredshift']) * xp.power(1+xp.array(events['redshift']),2.))
                prior = xp.exp(lnprior)

                inj_dict = {
                    'mass_1': xp.array(xp.array(events['mass1_source']) * (1 + xp.array(events['redshift']))),
                    'mass_2': xp.array(xp.array(events['mass2_source']) * (1 + xp.array(events['redshift']))),
                    'luminosity_distance': xp.array(events['luminosity_distance'])}
            
            else:
                raise ValueError('Catalog option not yet implemented. Please choose GWTC-4.0 or O3.')

            # If using the mass ratio, correct the prior with the Jacobian m2->q.
            if 'MassRatio' in pars['model-secondary']:
                inj_dict['mass_ratio'] = inj_dict.pop('mass_2') / inj_dict['mass_1']
                prior *= inj_dict['mass_1'] # |J_(m1,m2)->(m1,q)| = m1, with q = m2/m1.

        # Use simulated injections.
        else:

            print('\n\tUsing simulated injections to evaluate selection effects.')
            print('\n\t{}'.format(pars['injections-path']))

            try:
                if '.pickle' in pars['injections-path']:
                    with open(pars['injections-path'], 'rb') as f: data_inj = pickle.load(f)
                else:
                    raise ValueError('Only pickle files are currently supported for custom injections:\n{}'.format(pars['injections-path']))
            except:
                raise ValueError('Could not open the file containing the injections for selection effects. Please verify that the path is correct:\n{}'.format(pars['injections-path']))
            
            # This prior must be the one in detector frame for the variables (m1d,m2d,dL).
            # Whatever distribution and variables used to generate the injections, please make sure it follows such conventions.
            prior = xp.array(data_inj['prior'])
            obs_time = 1

            inj_dict = {
                'mass_1':              xp.array(data_inj['m1d']),
                'mass_2':              xp.array(data_inj['m2d']),
                'luminosity_distance': xp.array(data_inj['dL'])}
            
            if not pars['single-mass']:
                # If using the mass ratio, correct the prior with the Jacobian m2->q.
                if 'MassRatio' in pars['model-secondary']:
                    if not pars['inverse-mass-ratio']:
                        inj_dict['mass_ratio'] = inj_dict.pop('mass_2') / inj_dict['mass_1']
                        prior *= inj_dict['mass_1']                             # |J_(m1,m2)->(m1,q)| = m1, with q = m2/m1.
                    else:
                        inj_dict['mass_ratio'] = inj_dict['mass_1'] / inj_dict.pop('mass_2')
                        prior *= inj_dict['mass_1'] / inj_dict['mass_ratio']**2 # |J_(m1,m2)->(m1,q)| = m1/q^2, with q = m1/m2.
            else:
                # If only using one mass, remove the Jacobian contribution from the secondary.
                # This operation depends on the injection prior used to generate the injections.
                prior *= (1 + ref_cosmo.dl2z(inj_dict['luminosity_distance']))
                inj_dict.pop('mass_2')

        self.injections = icarogw.injections.injections(inj_dict, prior = prior, ntotal = pars['injections-number'], Tobs = obs_time)
        if not pars['catalog'] == 'GWTC-4.0':
            if   pars['selection-effects-cut'] == 'snr' : self.injections.update_cut(data_inj['snr'] >= pars['snr-cut' ])
            elif pars['selection-effects-cut'] == 'ifar': self.injections.update_cut(ifarmax         >= pars['ifar-cut'])
            else:
                raise ValueError('Unknown option to compute the selection effects cut.')
        else:
            self.injections.update_cut(selected_filt)

        print('\n\tUsing {} injections out of {} to compute selection effects.'.format(len(self.injections.injections_data['mass_1']), len(self.injections.injections_data_original['mass_1'])))

    def return_SelectionEffects(self):
        return self.injections



class Data:
       
    def __init__(self, pars):
        
        print('\n * Loading data.\n\n\t{}'.format(pars['data-path']))
    
        if not pars['true-data']:
            if   pars['PE-prior-distance'] == 'dL'   :     print('\n\tUsing a prior for PE samples uniform in luminosity distance.'               )
            elif pars['PE-prior-distance'] == 'dL3'  :     print('\n\tUsing a prior for PE samples uniform in comoving volume.'                   )
            elif pars['PE-prior-distance'] == 'UniformSourceFrame':     print('\n\tUsing a prior for PE samples uniform in source frame.'         )
            elif pars['PE-prior-distance'] == 'per-run':   print('\n\tUsing a prior for PE samples observing run-specific. See below for each event.')
            else:
                raise ValueError('Unknown option for PE sample prior distance.')
            if   pars['PE-prior-masses'  ] == 'm1-m2':     print('\n\tUsing a prior for PE samples uniform in component masses, (m1, m2).'        )
            elif pars['PE-prior-masses'  ] == 'Mc-q' :     print('\n\tUsing a prior for PE samples uniform in chirp mass and mass ratio, (Mc, q).')
            else:
                raise ValueError('Unknown option for PE sample prior masses.'  )
        else: print('\n\tIgnoring the PE priors as we just use the events true values.')

        if not pars['single-mass']:
            if 'MassRatio' in pars['model-secondary']:
                if not pars['inverse-mass-ratio']: print('\n\tUsing the mass ratio for the secondary, defined as q=m2/m1.')
                else:                              print('\n\tUsing the mass ratio for the secondary, defined as q=m1/m2.')
        else: print('\n\tUsing just the primary mass.')
        
        # Real GW data.
        if pars['real-data']:

            print('\n\tUsing IGWN catalogs of real GW events.')

            # FIXME: Implement option to select events based on FAR and/or SNR.
            # Need to add a new FAR/SNR key in the IGNW_events files.

            # The PE samples are publicly available here:
            # https://zenodo.org/records/6513631PE for O1, O2 and O3a events.
            # https://zenodo.org/records/8177023 for O3b events.
            # https://zenodo.org/records/17014085 for O4a events.

            from IGWN_pointers import O1_O2_BBHs_FAR_1, O3_BBHs_FAR_1, O4a_BBHs_FAR_1

            if   pars['catalog'] == 'GWTC-3':   catalog = O1_O2_BBHs_FAR_1 | O3_BBHs_FAR_1
            elif pars['catalog'] == 'GWTC-4.0': catalog = O1_O2_BBHs_FAR_1 | O3_BBHs_FAR_1 | O4a_BBHs_FAR_1
            elif pars['catalog'] == 'O3':       catalog = O3_BBHs_FAR_1
            elif pars['catalog'] == 'O4a':      catalog = O4a_BBHs_FAR_1
            else:
                raise ValueError('Unknown catalog option. Please choose from GWTC-3, GWTC-4.0, O3 or O4a.')

            print('')
            samps_dict = {}
            for ev in sorted(catalog.keys()):
                
                # Skip the events to be removed.
                if ev in pars['remove-events']: continue
                else:                           event_print = '\t{:<20}'.format(ev)

                tmp = h5py.File(catalog[ev]['PE'].replace("~IGNW_data_path", pars['data-path']))
                data_evs = tmp[catalog[ev]['PE_waveform']]['posterior_samples']

                pos_dict  = {
                    'mass_1'             : xp.array(data_evs['mass_1'][()]),
                    'mass_2'             : xp.array(data_evs['mass_2'][()]),
                    'luminosity_distance': xp.array(data_evs['luminosity_distance'][()])}

                # Account for PE priors. For O3 data, PE priors are uniform in component masses.
                # Luminosity distance.
                if   pars['PE-prior-distance'] == 'dL' : 
                    prior = xp.ones(len(pos_dict['luminosity_distance']))    # Set the prior to one.
                    event_print += " | dL prior: {:<30}".format('uniform in dL')

                elif pars['PE-prior-distance'] == 'dL3': 
                    prior = xp.power(   pos_dict['luminosity_distance'], 2.) # PE prior uniform in comoving volume: p(dL) \propto dL^2.
                    event_print += " | dL prior: {:<30}".format('uniform in detected volume')

                elif pars['PE-prior-distance'] == 'UniformSourceFrame':
                    prior_usf = bilby.gw.prior.UniformSourceFrame(
                        name='luminosity_distance', 
                        minimum=0.1, 
                        maximum=float(1.1*max(pos_dict['luminosity_distance'])), 
                        unit='Mpc'
                    )
                    prior = icarogw.cupy_pal.np2cp(prior_usf.prob(icarogw.cupy_pal.cp2np(pos_dict['luminosity_distance'])))
                    event_print += " | dL prior: {:<30}".format('uniform in source frame')

                elif pars['PE-prior-distance'] == 'per-run':
                    # If the filenames for O1-O3 events have 'nocosmo', they correspond to PE samples with dL priors \propto dL^2
                    if catalog[ev]['run'] in ['O1', 'O2', 'O3a', 'O3b'] and 'nocosmo' in catalog[ev]['PE']:
                        prior = xp.power(   pos_dict['luminosity_distance'], 2.)
                        event_print += " | dL prior: {:<30}".format('uniform in detected volume')
                    # If the filenames for O1-O3 events doesn't have 'nocosmo', they correspond to PE samples with UniformSourceFrame dL priors
                    elif catalog[ev]['run'] in ['O1', 'O2', 'O3a', 'O3b'] or catalog[ev]['run'] == 'O4a':
                        prior_usf = bilby.gw.prior.UniformSourceFrame(
                            name='luminosity_distance', 
                            minimum=0.1, 
                            maximum=float(1.1*max(pos_dict['luminosity_distance'])), 
                            unit='Mpc'
                        )
                        prior = icarogw.cupy_pal.np2cp(prior_usf.prob(icarogw.cupy_pal.cp2np(pos_dict['luminosity_distance'])))
                        event_print += " | dL prior: {:<30}".format('uniform in source frame')
                    else:
                        raise KeyError("Unknown run for event {} in IGWN_pointers dictionary.".format(catalog[ev]['run']))

                else:
                    raise ValueError("Unknown PE-prior-distance option. Please choose from 'dL', 'dL3', 'UniformSourceFrame', 'per-run'.")

                # Case of using mass ratio instead of the secondary mass.
                if 'MassRatio' in pars['model-secondary']:
                    pos_dict['mass_ratio'] = pos_dict.pop('mass_2') / pos_dict['mass_1']
                    prior *= pos_dict['mass_1'] # |J_(m1,m2)->(m1,q)| = m1, with q = m2/m1.
                
                samps_dict[ev] = icarogw.posterior_samples.posterior_samples(pos_dict, prior = prior)
                print(event_print)

        # Internal simulations.
        else:

            print('\n\tUsing a simulated catalog of GW events.')

            if '.pickle' in pars['data-path']:
                with open(pars['data-path'], 'rb') as f: data_evs = pickle.load(f)
            else:
                raise ValueError('Unknown format for the file containing the single events samples. Please make sure the file is correct:\n{}'.format(pars['data-path']))

            if pars['PE-prior-distance'] == 'per-run': raise ValueError("'per-run' PE-prior-distance option incompatible with simulated data.")

            samps_dict = {}
            for i in range(len(data_evs['m1d'])):
                idx = i

                if pars['true-data']:
                    pos_dict = {
                        'mass_1':              xp.array([data_evs['m1d'][i]]),
                        'mass_2':              xp.array([data_evs['m2d'][i]]),
                        'luminosity_distance': xp.array([data_evs['dL' ][i]])}
                else:
                    idx = list(data_evs['m1d'].keys())[i]
                    if idx in pars['remove-events']: continue
                    pos_dict = {
                        'mass_1':              xp.array(data_evs['m1d'][idx]),
                        'mass_2':              xp.array(data_evs['m2d'][idx]),
                        'luminosity_distance': xp.array(data_evs['dL' ][idx])}

                # Initialize the PE prior as flat for all variables. This is the case when only true values are used instead of the full PE.
                prior = xp.full(len(pos_dict['mass_1']), 1.)

                if 'MassRatio' in pars['model-secondary']:
                    if not pars['inverse-mass-ratio']: pos_dict['mass_ratio'] = pos_dict['mass_2'] / pos_dict['mass_1']
                    else:                              pos_dict['mass_ratio'] = pos_dict['mass_1'] / pos_dict['mass_2']

                # Account for PE prior.
                if not pars['true-data']:

                    # Luminosity distance. If the prior is uniform in dL, we leave it flat.
                    if pars['PE-prior-distance'] == 'dL3': 
                        prior *= data_evs['dL'][idx]**2 # PE prior uniform in comoving volume: p(dL) \propto dL^3.

                    if not pars['single-mass']:
                        chirp_mass = (pos_dict['mass_1'] * pos_dict['mass_2'])**(3./5.) / (pos_dict['mass_1'] + pos_dict['mass_2'])**(1./5.)
                        # Case of using component masses. If the prior is uniform in (m1,m2), we leave it flat.
                        if not 'MassRatio' in pars['model-secondary']:
                            if   pars['PE-prior-masses'] == 'Mc-q':
                                if not pars['inverse-mass-ratio']: prior *= chirp_mass / pos_dict['mass_1']**2 # |J_(Mc,q)->(m1,m2)| = Mc/m1^2, with q = m2/m1.
                                else:                              prior *= chirp_mass / pos_dict['mass_2']**2 # |J_(Mc,q)->(m1,m2)| = Mc/m2^2, with q = m1/m2.

                        else: # Case of using mass ratio instead of the secondary mass.
                            if not pars['inverse-mass-ratio']:
                                if   pars['PE-prior-masses'] == 'm1-m2': prior *= pos_dict['mass_1']              # |J_(m1,m2)->(m1,q)| = m1, with q = m2/m1.
                                elif pars['PE-prior-masses'] == 'Mc-q' : prior *= chirp_mass / pos_dict['mass_1'] # |J_(Mc,q)->(m1,q)| = Mc/m1, with q = m2/m1.
                            else:
                                if   pars['PE-prior-masses'] == 'm1-m2': prior *= pos_dict['mass_1'] / pos_dict['mass_ratio']**2 # |J_(m1,m2)->(m1,q)| = m1/q^2, with q = m1/m2.
                                elif pars['PE-prior-masses'] == 'Mc-q' : prior *= chirp_mass / pos_dict['mass_1']                # |J_(Mc,q)->(m1,q)| = Mc/m1, with q = m1/m2.

                samps_dict['{}'.format(i)] = icarogw.posterior_samples.posterior_samples(pos_dict, prior = prior)
        
        self.data = icarogw.posterior_samples.posterior_samples_catalog(samps_dict)
        print('\n\tUsing a population of {} events.'.format(self.data.n_ev))

    def return_Data(self):
        return self.data



class LikelihoodPrior:

    def __init__(self, pars, data, injections, wrapper):
          
        self.pars       = pars
        self.data       = data
        self.injections = injections
        self.wrapper    = wrapper

    def Likelihood(self, data, injections, wrapper):

        if self.pars['loglike-var'] == 0: self.pars['loglike-var'] = None
        
        if not self.pars['ignore-selection-effects']:
            res = icarogw.likelihood.hierarchical_likelihood(
                            data, injections, wrapper,
                            nparallel               = self.pars['nparallel'],
                            neffPE                  = self.pars['neffPE'],
                            neffINJ                 = self.pars['neffINJ'],
                            likelihood_variance_thr = self.pars['loglike-var'])
        else: # Use the likeliood without selection effects.
            res = icarogw.likelihood.hierarchical_likelihood_no_selection_effects(
                            data, wrapper,
                            nparallel               = self.pars['nparallel'],
                            neffPE                  = self.pars['neffPE'])
        return res

    def Prior(self, pars, w):

        def initialise_prior(dict_in, dict_out, w):

            available_bilby_priors = [
                'Uniform',
                'LogUniform',
            ]

            # Precompute spline coefficient names safely
            use_dirichlet = ('Spline' in pars['model-primary']) and pars['dirichlet-prior']
            if use_dirichlet:
                spline_coeffs = [f'c{i}' for i in range(1, pars['splines-number']+1)]
              
            for par in w.population_parameters:

                # Dirichlet case via Gamma priors
                if use_dirichlet and par in spline_coeffs:
                    dict_out[par] = bilby.core.prior.Gamma(1., 1., name = par)

                # Standard Uniform / LogUniform priors
                elif isinstance(dict_in[par], list) and len(dict_in[par]) == 2:
                    dict_out[par] = bilby.core.prior.Uniform(dict_in[par][0], dict_in[par][1])

                elif isinstance(dict_in[par], list) and len(dict_in[par]) > 2:
                    if dict_in[par][2] in available_bilby_priors:
                        bilby_prior_class = getattr(bilby.core.prior, dict_in[par][2])
                        dict_out[par] = bilby_prior_class(dict_in[par][0], dict_in[par][1])
                    else:
                        raise KeyError("Unknown bilby prior. Available (in this pipeline):\n\t" + '\n\t'.join(available_bilby_priors))

                # Fixed value
                elif isinstance(dict_in[par], float): dict_out[par] = dict_in[par]

                else:
                    raise ValueError("Unknown type for prior on {}. Please provide either a fixed value, or a 2-list [min, max], or a 3-list [min, max, type]".format(dict_in[par]))
            
            print('\n * Using the following priors.\n')
            if use_dirichlet:
                print_dictionary({key: dict_in[key] for key in dict_out.keys()-set(spline_coeffs)})
                print('\n\tUsing a Dirichlet prior for the spline coefficients: {}.'.format(spline_coeffs))
            else:
                print_dictionary({key: dict_in[key] for key in dict_out.keys()})

            print('\n\tWith constraints:')
            # Miscellaneous built-in additional constraints for some models
            constraints_dict = {
                'MD_redundancy': {
                    'pars':        ['gamma', 'kappa'], 
                    'func':        (lambda x, y: x + y),
                    'const_bilby': bilby.core.prior.Constraint(minimum = 0, 
                                                               maximum = xp.inf),
                    'print':       "\t[ gamma + kappa > 0 ]",
                },
                'w0wa_earlyMDera': {
                    'pars':        ['w0', 'wa'], 
                    'func':        (lambda x, y: x + y),
                    'const_bilby': bilby.core.prior.Constraint(minimum = - xp.inf, 
                                                               maximum = 0),
                    'print':       "\t[ w0 + wa < 0 ]",
                },
                'MLTP_peak_ordering': {
                    'pars':        ['mu_g_high', 'mu_g_low'], 
                    'func':        (lambda x, y: x - y),
                    'const_bilby': bilby.core.prior.Constraint(minimum = 0., 
                                                               maximum = xp.inf),
                    'print':       "\t[ mu_g_high > mu_g_low ]",
                },
                'PL2G_peak_ordering': {
                    'pars':        ['mu_b_z0', 'mu_a_z0'], 
                    'func':        (lambda x, y: x - y),
                    'const_bilby': bilby.core.prior.Constraint(minimum = 0., 
                                                               maximum = xp.inf),
                    'print':       "\t[ mu_b_z0 > mu_a_z0 ]",
                },
                'nPL_peak_ordering': {
                    'pars':        [f'mmin_{c}' for c in "abcdefghij"], 
                    'func':        lambda *mmins: minimum.reduce([ (mmins[i+1] - mmins[i]) for i in range(len(mmins)-1)]),
                    'const_bilby': bilby.core.prior.Constraint(minimum = 0., 
                                                               maximum = xp.inf),
                    'print':       "\t[ mmin_{p} < mmin_{p+1} for all peaks p ] (PL peak ordering).",
                },
                'nPL_minmax_ordering': {
                    'pars':        [f'{p}_{c}' for c in "abcdefghij" for p in ['mmin', 'mmax']], 
                    'func':        lambda *pars: minimum.reduce([ (pars[2*i+1] - pars[2*i]) for i in range(len(pars)//2)]),
                    'const_bilby': bilby.core.prior.Constraint(minimum = 1., 
                                                               maximum = xp.inf),
                    'print':       "\t[ mmin_{p} < mmax_{p} for peaks p ] (PL minmax ordering).",
                },
            }

            # MD redundancy constraint
            if not (pars['model-rate'] == 'MadauDickinson' and pars['constraint_MD_redundancy']):
                constraints_dict.pop('MD_redundancy')
            elif pars['model-rate'] != 'MadauDickinson' and pars['constraint_MD_redundancy']:
                raise ValueError("MD redundancy constraint is only available when using MadauDickinson rate evolution parametrization.")
            else:
                pass
            # w0wa early matter domination era constraint
            if not (pars['model-cosmology'] == 'Flatw0waCDM' and pars['constraint_w0wa_earlyMDera']):
                constraints_dict.pop('w0wa_earlyMDera')
            elif pars['model-cosmology'] != 'Flatw0waCDM' and pars['constraint_w0wa_earlyMDera']:
                raise ValueError("w0wa early MD era constraint is only available when using Flatw0waCDM cosmological model.")
            else:
                pass
            # peak ordering constraints
            if pars['constraint_peak_ordering']:
                # MLTP
                if not pars['model-primary'] == 'PowerLaw-Gaussian-Gaussian': 
                    constraints_dict.pop('MLTP_peak_ordering')
                else:
                    pass
                # PL2G
                if pars['model-primary'] == 'PowerLaw-GaussianRedshiftLinear-GaussianRedshiftLinear':
                    # must not include redshift evolution
                    if (
                            (pars['redshift-mixture']) or
                            (not (
                                type(dict_in['mu_a_z1']) == float and
                                type(dict_in['mu_b_z1']) == float and
                                type(dict_in['sigma_a_z1']) == float and
                                type(dict_in['sigma_b_z1']) == float
                            ))
                    ): 
                        raise ValueError("peak ordering constraint for 'PowerLaw-GaussianRedshiftLinear-GaussianRedshiftLinear' mass model is only available with no redshift evolution.")
                    else:
                        pass
                else:
                    constraints_dict.pop('PL2G_peak_ordering')

                # 2PL, 3PL, 4PL (also implements minmax ordering for all peaks)
                if pars['model-primary'] == 'PowerLaw-PowerLaw': 
                    constraints_dict['nPL_peak_ordering']['pars']   = constraints_dict['nPL_peak_ordering']['pars'][:2]
                    constraints_dict['nPL_minmax_ordering']['pars'] = constraints_dict['nPL_minmax_ordering']['pars'][:2*2]
                elif pars['model-primary'] == 'PowerLaw-PowerLaw-PowerLaw': 
                    constraints_dict['nPL_peak_ordering']['pars']   = constraints_dict['nPL_peak_ordering']['pars'][:3]
                    constraints_dict['nPL_minmax_ordering']['pars'] = constraints_dict['nPL_minmax_ordering']['pars'][:3*2]
                elif pars['model-primary'] == 'PowerLaw-PowerLaw-PowerLaw-PowerLaw': 
                    constraints_dict['nPL_peak_ordering']['pars']   = constraints_dict['nPL_peak_ordering']['pars'][:4]
                    constraints_dict['nPL_minmax_ordering']['pars'] = constraints_dict['nPL_minmax_ordering']['pars'][:4*2]
                else: 
                    constraints_dict.pop('nPL_peak_ordering')
                    constraints_dict.pop('nPL_minmax_ordering')
            else:
                constraints_dict.pop('MLTP_peak_ordering')
                constraints_dict.pop('PL2G_peak_ordering')
                constraints_dict.pop('nPL_peak_ordering')
                constraints_dict.pop('nPL_minmax_ordering')

            # implementing the conversion function based on all the constraints that we kept
            def constraints_conversion_function(params):
                converted_params = params.copy()
                for const in constraints_dict:
                    converted_params[const] = constraints_dict[const]['func'](
                        *[params[par] for par in constraints_dict[const]['pars']]
                    )
                return converted_params

            # Adding the constraitns as bilby objects in the prior dictionary
            for const in constraints_dict:
                dict_out[const] = constraints_dict[const]['const_bilby']
                print(constraints_dict[const]['print'])
            
            dict_out = bilby.core.prior.PriorDict(dict_out, conversion_function = constraints_conversion_function)

            return dict_out

        prior = bilby.core.prior.PriorDict()
        prior = initialise_prior(pars['all-priors'], prior, w)

        return prior
        
    def return_LikelihoodPrior(self):

        self.likelihood = self.Likelihood(self.data, self.injections, self.wrapper)
        self.prior      = self.Prior(self.pars, self.wrapper)

        return self.likelihood, self.prior



def main():

    mp.set_start_method("fork", force=True)

    # ------------------------------------------------ #
    # Read config file and initialise input parameters #
    # ------------------------------------------------ #

    parser = OptionParser(options.usage)
    parser.add_option(      '--config-file', type='string', metavar = 'config_file', default = None)
    parser.add_option('-n', '--n-processes', type='int',    metavar = 'n_processes', default = -1, help="Set the number of processes for parallelized injections generation from command line, if this should match some external structure (e.g. number of CPUs allocated to the simulation on a computing cluster job.)")
    (opts, _) = parser.parse_args()

    config_file = opts.config_file
    if not config_file: parser.error('Please specify a config file.\n')
    if not os.path.exists(config_file): parser.error('Config file {} not found.\n'.format(config_file))

    Config = configparser.ConfigParser()
    Config.read(config_file)

    # Initialise input parameters dictionary.
    input_pars = options.InitialiseOptions(Config)

    # Set output directory.
    if not os.path.exists(input_pars['output']): os.makedirs(input_pars['output'])

    # Copy config file to output.
    try:    shutil.copyfile(config_file, os.path.join(input_pars['output'], os.path.basename(os.path.normpath(config_file))))
    except: pass # Config file already copied.

    # Set the number of parallel processes according to command line if provided to match hardware structure
    if opts.n_processes > 0:
        print(f"\n * Number of processes set via command-line option: n_processes = {opts.n_processes} \n")
        input_pars['npool'] = opts.n_processes

    # Deviate stdout and stderr to file.
    if not input_pars['screen-output']:
        print("\n * screen-output is False -> I will deviate stdout and stderr to the result folder.\n")
        sys.stdout = open(os.path.join(input_pars['output'], 'stdout_icarogw.txt'), 'w')
        sys.stderr = open(os.path.join(input_pars['output'], 'stderr_icarogw.txt'), 'w')
    else: pass
    print('\n\n ===== Starting  i c a r o g w  runner ===== \n')

    if icarogw.cupy_pal.is_there_cupy(): print(" * I will run with cupy on GPU.\n")

    # Print run parameters.
    print(' * I will be running with the following parameters.\n')
    print_dictionary(input_pars)

    # ------------------------------------------------------------------------------ #
    # Initialise the ICAROGW run: model, data, selection effects, likelihood, priors #
    # ------------------------------------------------------------------------------ #

    # Initialise the model wrappers.
    tmp = Wrappers(input_pars)
    m1w, m2w, rw, cw, ref_cosmo = tmp.return_Wrappers()

    tmp = Rate(input_pars, m1w, m2w, rw, cw)
    wrapper = tmp.return_Rate()

    # Read injections for selection effects.
    if not input_pars['ignore-selection-effects']:
        tmp = SelectionEffects(input_pars, ref_cosmo)
        injections = tmp.return_SelectionEffects()
    else:
        injections = None

    # Read events data.
    tmp = Data(input_pars)
    data = tmp.return_Data()

    # Initialise hierarchical likelihood and set the priors.
    tmp = LikelihoodPrior(input_pars, data, injections, wrapper)
    likelihood, prior = tmp.return_LikelihoodPrior()

    # Plot weighted injections
    # print("\n * Plotting weighted injections")
    # try: 
    #     postprocessing.plot_weighted_injections(input_pars, injections=injections, rate=wrapper, data=data)
    #     print("\t...done !")
    # except: 
    #     pass
    #     print("\t...failed. Carry on...\n")

    if not input_pars['ignore-selection-effects']:
        # Control the effective number of injections on the injected model.
        check_effective_number_injections(input_pars, likelihood, data.n_ev)

    # ----------------------------------------------- #
    # Start the sampler and run hierarchical analysis #
    # ----------------------------------------------- #

    # Extracting the input parameters specific to each class of sampler 
    print('\n * Running hierarchical analysis with this settings.\n')
    if   input_pars['sampler'] == 'dynesty' or input_pars['sampler'] == 'nessai':
        sampler_pars = {key: input_pars[key] for key in ['sampler', 'nlive', 'naccept', 'print-method', 'sample']}
        if icarogw.cupy_pal.is_there_cupy():
            sampler_pars['npool'] = None
            print("\tRunning sampler on GPU, enforcing npool==None (no multiprocessing).")
        else: 
            sampler_pars['npool'] = input_pars['npool']
        print_dictionary(sampler_pars)

    elif input_pars['sampler'] == 'ptemcee':
        sampler_pars = {key: input_pars[key] for key in ['sampler', 'nwalkers', 'ntemps', 'threads', 'print-method']}
        print_dictionary(sampler_pars)

    elif input_pars['sampler'] == 'emcee': 
        sampler_pars = {key: input_pars[key] for key in ['sampler', 'nwalkers', 'nsteps', 'npool']}
        print_dictionary(sampler_pars)

    else:
        raise ValueError('Sampler not available.')

    if input_pars['sampler'] == 'nessai':
        if not input_pars['nessai-plot']: sampler_pars.update(dict(nessai_plot = False))

    # Start Bilby sampler.
    print('\n * Starting the sampler.\n')
    hierarchical = bilby.run_sampler(
            likelihood, prior,
            outdir       = os.path.join(input_pars['output'], 'sampler'),
            **sampler_pars,
    )
    if input_pars['true-values'] == {}: hierarchical.plot_corner()
    else:                               hierarchical.plot_corner(truth = {key: input_pars['true-values'][key] for key in hierarchical.search_parameter_keys})
    
    # Get the samples.
    samp_path = os.path.join(input_pars['output'], 'sampler', 'label_result.json')
    with open(samp_path) as f:
        tmp = json.load(f)
        df  = pd.DataFrame(tmp['posterior']['content'])
        priors_dict = tmp['priors']

    # Save the evidence.
    with open('{}/log_evidence.txt'.format(input_pars['output']), 'w') as f:
        f.write('{}\n'.format('# log_Z_base_e\tlog_Z_err\tmax_log_L'))
        log_evidence_err = round(tmp['log_evidence_err'], 2)
        if xp.isnan(tmp['log_evidence_err']): log_evidence_err = 0.1
        f.write('{}\t{}\t\t{}'.format(round(tmp['log_evidence'], 2), log_evidence_err, round(max(df['log_likelihood']), 2)))

    # Control the effective number of injections on the maximum likelihood model.
    print('\n * Computing effective number of injections.')
    maxL_index  = int(xp.argmax(xp.array(df['log_likelihood'])))
    maxL_values = {key: df[key][maxL_index] for key in wrapper.population_parameters}
    if not input_pars['ignore-selection-effects']:
        # Control the effective number of injections on the maximum likelihood model.
        print('\n * Computing effective number of injections.')
        check_effective_number_injections(input_pars, likelihood, data.n_ev, maxL_values = maxL_values)
    print('\n * Maximum likelihood values.\n')
    print_dictionary(maxL_values)

    # ----------------------------------- #
    # Plots production and postprocessing #
    # ----------------------------------- #

    print('\n * Producing plots.')
    input_pars['output-plots']  = os.path.join(input_pars['output'], 'plots' )
    if not os.path.exists(input_pars['output-plots']):  os.makedirs(input_pars['output-plots'] )

    tmp = postprocessing.Plots(input_pars, df, m1w, m2w, rw, cw, ref_cosmo, wrapper, priors_dict, injections)
    tmp.ProducePlots()

    # Save curves of the reconstructed distributions.
    curves = tmp.return_curves()
    with open(os.path.join(input_pars['output-plots'], 'curves.pkl'), 'wb') as f:
        pickle.dump(curves, f)

    print('\n * Finished.\n')
