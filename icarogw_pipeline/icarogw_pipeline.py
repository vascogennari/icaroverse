import os, sys, configparser, shutil
from optparse import OptionParser
from options import usage

import pickle, json, h5py
import pandas as pd, numpy as np
import icarogw, bilby


def config_generator(base_pars, pars, key, item):

    config_file = os.path.join(base_pars['output-directory'], 'name.ini')
    shutil.copyfile(base_pars['default-config'], config_file)
    config = configparser.ConfigParser()
    config.read(config_file)

    # Input
    #config.set('input', 'output',          '{}'.format(pars['output']         ))
    #config.set('input', 'injections-path', '{}'.format(pars['injections-path']))
    #config.set('input', 'data-path',       '{}'.format(pars['data-path']      ))

    config.set('input', '{}'.format(key), '{}'.format(item))

    with open(config_file, 'w') as f: config.write(f)


def Selection():

    def Models():

        tmp = [
            'PowerLaw-GaussianRedshiftLinear'
        ]
        return tmp
    
    def Data():

        tmp = [
        ]
        return tmp

    
    def Injections():

        tmp = [
        ]
        return tmp
    
    selection_dictionary = {
        'model-primary'   : Models(),
        'data-path'       : Data(),
        'injections-path' : Injections(),
    }

    return selection_dictionary


def main():

    print('\n\n Starting  i c a r o g w  pipeline \n')

    '''
        Generates multiple config files
        Generates and launches condor files
    '''

    input_pars  = {
        'output-directory' : 'configs/TEST',
        'default-config'   : 'configs/config_default.ini',
    }

    config_pars = {
        'output'           : '/Users/vgennari/Desktop/INJ-evolving_REC-evolving',
        'injections-path'  : '/Users/vgennari/Documents/work/code/python/icarogw/data/simulations/injections_selection_effects/inj_PowerLawPeak_N10000000_SNR12_fGW15/inj_PowerLawPeak_N10000000_SNR12_fGW15.pickle',
        'data-path'        : '/Users/vgennari/Documents/work/code/python/icarogw/data/simulations/simulated_population/pop-252600_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio_PowerLaw_evolving/events_detector_dict_pop-252600_PowerLawRedshiftLinear-GaussianRedshiftLinear_MassRatio_PowerLaw_evolving.pickle',
    }

    if not os.path.exists(input_pars['output-directory']): os.makedirs(input_pars['output-directory'])

    select_dict = Selection()

    for key in select_dict.keys():
        for item in select_dict[key]:
            config_generator(input_pars, config_pars, key, item)

    print('\n * Finished.\n')

if __name__=='__main__':
    main()
