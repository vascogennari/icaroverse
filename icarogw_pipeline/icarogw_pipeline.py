import os, configparser, shutil


def config_initialiser(base_pars, rate, secondary, primary, extra_text):

    name = '{}_{}_{}_{}_{}'.format(base_pars['analysis_name'], primary, secondary, rate, extra_text)
    print('\t{}'.format(name))
    config_path = os.path.join(base_pars['configs-directory'], 'config_{}.ini'.format(name))
    shutil.copyfile(base_pars['default-config'], config_path)
    return configparser.ConfigParser(), config_path, name


def config_writer(config, config_path, pars, inj, data, rate, secondary, primary, name):

    config.read(config_path)

    # input
    config.set('input', 'output', '{}/{}'.format(pars['output'], name))
    config.set('input', 'injections-path', inj)
    config.set('input', 'data-path',       data)

    for key in ['simluation', 'injections-number', 'selection-effects-cut', 'snr-cut']:
        if key in pars.keys(): config.set('input',   key, str(pars[key]))

    # model
    config.set('model', 'model-primary',   primary)
    config.set('model', 'model-secondary', secondary)
    config.set('model', 'model-rate',      rate)

    for key in ['redshift-transition', 'low-smoothing', 'scale-free', 'priors']:
        if key in pars.keys(): config.set('model',   key, str(pars[key]))

    # sampler
    for key in ['nparallel', 'neffPE', 'npool', 'nlive']:
        if key in pars.keys(): config.set('sampler', key, str(pars[key]))

    # plots
    for key in ['true-values', 'selection-effects']:
        if key in pars.keys(): config.set('plots',   key, str(pars[key]))

    with open(config_path, 'w') as f: config.write(f)


def filter_selection(dictionary):

    return [key for key, value in dictionary.items() if value == True]


def Selection():

    def Model_Primary():
        tmp = {
            'PowerLaw-GaussianRedshiftLinear':  1,
            'GaussianRedshiftLinear-order-1' :  0,
        }
        return filter_selection(tmp)
    
    def Model_Secondary():
        tmp = {
            'Mass2-PowerLaw'                 :  0,
            'MassRatio-Gaussian'             :  1,
            'MassRatio-PowerLaw'             :  0,
        }
        return filter_selection(tmp)

    def Model_Rate():
        tmp = {
            'PowerLaw'                       :  1,
            'MadauDickinson'                 :  0,
        }
        return filter_selection(tmp)
        
    def Data():
        tmp = {
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_1/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_1.pickle'  : 'inj-NNN-1_rec-NNE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_2/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_2.pickle'  : 'inj-NNN-2_rec-NNE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_3/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_3.pickle'  : 'inj-NNN-3_rec-NNE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_4/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_4.pickle'  : 'inj-NNN-4_rec-NNE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_5/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_5.pickle'  : 'inj-NNN-5_rec-NNE_low-SNR',
            
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_1/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_1.pickle': 'inj-NNN-1_rec-NNE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_2/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_2.pickle': 'inj-NNN-2_rec-NNE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_3/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_3.pickle': 'inj-NNN-3_rec-NNE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_4/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_4.pickle': 'inj-NNN-4_rec-NNE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_5/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_5.pickle': 'inj-NNN-5_rec-NNE_high-SNR',

            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_1/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_1.pickle'  : 'inj-NEN-1_rec-NNE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_2/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_2.pickle'  : 'inj-NEN-2_rec-NNE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_3/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_3.pickle'  : 'inj-NEN-3_rec-NNE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_4/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_4.pickle'  : 'inj-NEN-4_rec-NNE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_5/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_5.pickle'  : 'inj-NEN-5_rec-NNE_low-SNR',

            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_1/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_1.pickle': 'inj-NEN-1_rec-NNE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_2/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_2.pickle': 'inj-NEN-2_rec-NNE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_3/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_3.pickle': 'inj-NEN-3_rec-NNE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_4/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_4.pickle': 'inj-NEN-4_rec-NNE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_5/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_5.pickle': 'inj-NEN-5_rec-NNE_high-SNR',

            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_1/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_1.pickle'  : 'inj-NNN-1_rec-NEE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_2/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_2.pickle'  : 'inj-NNN-2_rec-NEE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_3/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_3.pickle'  : 'inj-NNN-3_rec-NEE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_4/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_4.pickle'  : 'inj-NNN-4_rec-NEE_low-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_5/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_5.pickle'  : 'inj-NNN-5_rec-NEE_low-SNR',
            
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_1/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_1.pickle': 'inj-NNN-1_rec-NEE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_2/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_2.pickle': 'inj-NNN-2_rec-NEE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_3/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_3.pickle': 'inj-NNN-3_rec-NEE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_4/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_4.pickle': 'inj-NNN-4_rec-NEE_high-SNR',
            # '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_5/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NNN_5.pickle': 'inj-NNN-5_rec-NEE_high-SNR',

            '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_1/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_1.pickle'  : 'inj-NEN-1_rec-NEE_low-SNR',
            '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_2/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_2.pickle'  : 'inj-NEN-2_rec-NEE_low-SNR',
            '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_3/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_3.pickle'  : 'inj-NEN-3_rec-NEE_low-SNR',
            '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_4/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_4.pickle'  : 'inj-NEN-4_rec-NEE_low-SNR',
            '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_5/events_detector_dict_pop-18000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_5.pickle'  : 'inj-NEN-5_rec-NEE_low-SNR',

            '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_1/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_1.pickle': 'inj-NEN-1_rec-NEE_high-SNR',
            '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_2/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_2.pickle': 'inj-NEN-2_rec-NEE_high-SNR',
            '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_3/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_3.pickle': 'inj-NEN-3_rec-NEE_high-SNR',
            '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_4/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_4.pickle': 'inj-NEN-4_rec-NEE_high-SNR',
            '/home/vasco.gennari/icarogw/data/simulations/simulated_population/simulation_evolving_PROD1/pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_5/events_detector_dict_pop-350000_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio-Gaussian_PowerLaw_NEN_5.pickle': 'inj-NEN-5_rec-NEE_high-SNR',
        }
        return tmp

    
    def Injections():
        tmp = [
            '/home/vasco.gennari/icarogw/data/simulations/injections_selection_effects/inj_PowerLawPeak_N100000_SNR12_fGW15/inj_PowerLawPeak_N100000_SNR12_fGW15.pickle',
        ]
        return tmp
    
    selection_dictionary = {
        'model-primary'   : Model_Primary(),
        'model-secondary' : Model_Secondary(),
        'model-rate'      : Model_Rate(),
        'data-path'       : Data(),
        'injections-path' : Injections(),
    }

    return selection_dictionary


def main():

    print('\n\n Starting  i c a r o g w  pipeline')

    '''
        Generates multiple config files.
        Generates and launches condor files.
    '''

    config_pars  = {
        'configs-directory'     : '/Users/vgennari/Documents/work/code/python/icarogw_pipeline/config_files/simulation_evolving_population_PROD1',
        'default-config'        : '/Users/vgennari/Documents/work/code/python/icarogw_pipeline/config_files/config_default.ini',
        'analysis_name'         : 'simluation-evolution-PROD1'
    }

    input_pars   = {
        # input
        'output'                : '/home/vasco.gennari/public_html/icarogw/simulations/evolving_population/PROD1',
        'simulation'            : 1,
        'injections-number'     : 4430000000,
        'selection-effects-cut' : 'snr',
        'snr-cut'               : 12,

        # model
        'redshift-transition'   : 'linear',
        'low-smoothing'         : 1,
        'scale-free'            : 1,
        # 'priors'                : {'mu_z0': [20., 60.], 'mu_z1': 0., 'sigma_z0': [1. , 30. ], 'sigma_z1': 0., 'gamma': [-20., 20.], 'mu_q': [0.1, 1.0], 'sigma_q': [0.01, 0.9], 'mix_z1': 0.9},      # rec-NNE
        'priors'                : {'mu_z0': [20., 60.], 'mu_z1': [-80., 80.], 'sigma_z0': [1. , 30. ], 'sigma_z1': 0., 'gamma': [-20., 20.], 'mu_q': [0.1, 1.0], 'sigma_q': [0.01, 0.9], 'mix_z1': 0.9},      # rec-NEE

        # sampler
        'nparallel'             : 1,
        'neffPE'                : -1,
        'npool'                 : 20,
        'nlive'                 : 1000,

        # plots
        #  'true-values'           : {'H0': 67.7, 'Om0': 0.308, 'alpha': 3.8, 'mmin': 7.0, 'mmax': 150.0, 'mu_z0': 35.0, 'mu_z1':  0.0, 'sigma_z0': 6.0, 'sigma_z1': 0.0, 'mix_z0': 0.9, 'mix_z1': 0.9, 'delta_m': 5.0, 'mu_q': 0.8, 'sigma_q': 0.15, 'gamma': 0.0},      # inj-NNN
        'true-values'           : {'H0': 67.7, 'Om0': 0.308, 'alpha': 3.8, 'mmin': 7.0, 'mmax': 150.0, 'mu_z0': 35.0, 'mu_z1': 30.0, 'sigma_z0': 6.0, 'sigma_z1': 0.0, 'mix_z0': 0.9, 'mix_z1': 0.9, 'delta_m': 5.0, 'mu_q': 0.8, 'sigma_q': 0.15, 'gamma': 0.0},      # inj-NEN
        # 'true-values'           : {'H0': 67.7, 'Om0': 0.308, 'alpha': 3.8, 'mmin': 7.0, 'mmax': 150.0, 'mu_z0': 35.0, 'mu_z1':  0.0, 'sigma_z0': 6.0, 'sigma_z1': 0.0, 'mix_z0': 0.9, 'mix_z1': 0.9, 'delta_m': 5.0, 'mu_q': 0.8, 'sigma_q': 0.15, 'gamma': 3.0},      # inj-NNE
        # 'true-values'           : {'H0': 67.7, 'Om0': 0.308, 'alpha': 3.8, 'mmin': 7.0, 'mmax': 150.0, 'mu_z0': 35.0, 'mu_z1': 30.0, 'sigma_z0': 6.0, 'sigma_z1': 0.0, 'mix_z0': 0.9, 'mix_z1': 0.9, 'delta_m': 5.0, 'mu_q': 0.8, 'sigma_q': 0.15, 'gamma': 3.0},      # inj-NEE
        'selection-effects'     : 1,
    }
    # Remove from the input parameters those set to False.
    input_pars = {key: input_pars[key] for key, value in input_pars.items() if not value == False}

    if not os.path.exists(config_pars['configs-directory']): os.makedirs(config_pars['configs-directory'])
    print('\n * Config files are saved here.\n\n\t{}'.format(config_pars['configs-directory']))
    print('\n * Generating the following config files.\n')
    select_dict = Selection()
    count = 0

    # Loop on all the selected config files.
    for inj in select_dict['injections-path']:
        for data in select_dict['data-path'].keys():
            for rate in select_dict['model-rate']:
                for secondary in select_dict['model-secondary']:
                    for primary in select_dict['model-primary']:
                        # Initialise the config file.
                        config, config_path, name = config_initialiser(config_pars, rate, secondary, primary, select_dict['data-path'][data])
                        # Write on the config file.
                        config_writer(config, config_path, input_pars, inj, data, rate, secondary, primary, name)
                        count += 1

    print('\n * The config files created is {}'.format(count))
    print('\n * Finished.\n\n')

if __name__=='__main__':
    main()
