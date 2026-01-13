import os, sys, configparser
from optparse import OptionParser


template = """#!/bin/sh

#SBATCH --job-name={name}
#SBATCH --output={slurm_files_dir_path}/output_{name}.out
#SBATCH --error={slurm_files_dir_path}/error_{name}.err
#SBATCH --nodes={nodes}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory}G
#SBATCH --partition=htc
#SBATCH --time={time}
#SBATCH --account=virgo
#SBATCH --licenses=sps
{email_option}

module load conda
conda activate {conda_env}
{executable} {script} --config-file {config} -n $SLURM_CPUS_PER_TASK
"""

email_option_template = """
#SBATCH --mail-user={user_mail}
#SBATCH --mail-type=ALL
"""

def activate_slurm_submit(pars):

    sys.stderr.write('generating {}\n'.format(pars['submission_filepath']))

    if not os.path.exists(os.path.dirname(pars['slurm_files_dir_path'])): os.makedirs(os.path.dirname(pars['slurm_files_dir_path']))
    if not os.path.exists(pars['slurm_files_dir_path']):                  os.makedirs(pars['slurm_files_dir_path'])

    with open(pars['submission_filepath'], 'w') as f:
        if pars['user_mail'] != '': email_option = email_option_template.format(user_mail=pars['user_mail'])
        else: email_option = ''
        submission_command = template.format(
            name                = pars['job_name'],
            slurm_files_dir_path = pars['slurm_files_dir_path'],
            nodes               = pars['slurm_nodes'],
            cpus                = pars['slurm_cpus'],
            memory              = pars['slurm_memory'],
            time                = '{}-{}:{}:00'.format(pars['slurm_time']['days'], pars['slurm_time']['hours'], pars['slurm_time']['minutes']),
            email_option        = email_option,
            conda_env           = pars['conda_env'],
            executable          = pars['slurm_python_path'],
            script              = pars['slurm_executable_file'],
            config              = pars['config_filepath']
        )

        f.write(submission_command)
    
    sys.stderr.write('submitting {}\n\n'.format(pars['submission_filepath']))
    if not pars['test']: os.system('sbatch {}'.format(pars['submission_filepath']))
    else: pass

# ---------------------------------------------------------------------- #

def main():

    parser = OptionParser()
    parser.add_option('-p', '--population', action="store_true", dest='population', default=False, help="Flag to run the slurm pipeline on a population directory rather than a config files directory.")
    parser.add_option('-t', '--test',       action="store_true", dest='test',       default=False, help="Flag to test the slurm pipeline without actually launching slurm jobs.")
    opts, _ = parser.parse_args()

    pars['test'] = opts.test

    if opts.population:

        pars['population_PE_dir_path'] =  os.path.join(pars['population_dir_path'], "parameter_estimation")
        if not os.path.exists(pars['population_PE_dir_path']):
            raise FileNotFoundError("No `parameter_estimation` directory found in the population directory. Please make sure the `population_directory` path is correct and/or the corresponding population has been processed with `generate_configs.py`")

        events_dir_list = sorted([ev_dir for ev_dir in os.listdir(pars['population_PE_dir_path']) if 'event' in ev_dir])

        print('')
        for event_dir_name in events_dir_list:
            pars['event_dir_path'] = os.path.join(pars['population_PE_dir_path'], event_dir_name)
            config_filename = [fname for fname in os.listdir(pars['event_dir_path']) if 'config' in fname][0]

            # Entries needed for activate_slurm_submit()
            pars['config_filepath'] = os.path.join(pars['event_dir_path'], config_filename)
            pars['job_name'] = '_'.join([os.path.basename(pars['population_dir_path']), os.path.basename(pars['event_dir_path'])])
            pars['slurm_files_dir_path'] = pars['event_dir_path']
            pars['submission_filepath'] = os.path.join(pars['slurm_files_dir_path'], f"submit_{pars['job_name']}.sh")
            activate_slurm_submit(pars)

        print('\nThe PE config files in {configs_path} are running in {n_jobs} detached slurm jobs, within {conda_env} conda environment. Good luck!\n'.format(configs_path = pars['population_PE_dir_path'], conda_env=pars['conda_env'], n_jobs=len(events_dir_list)))

    else:

        config_filenames = [name for name in os.listdir(pars['config_files_dir_path']) if '.ini' in name]

        for config_filename in config_filenames:
            # Entries needed for activate_slurm_submit()
            pars['config_filepath'] = os.path.join(pars['config_files_dir_path'], config_filename)
            
            Config = configparser.ConfigParser()
            Config.read(pars['config_filepath'])
            try: pars['output_directory'] = Config.get('input', 'output')
            except: raise ValueError("Please make sure the config file {} has an `output` directory entry.")

            job_name = config_filename.strip(".ini").strip("config_")
            pars['job_name'] = job_name
            pars['slurm_files_dir_path'] = os.path.join(pars['output_directory'], "slurm_files")
            submission_filename = f"submit_{pars['job_name']}.sh"
            pars['submission_filepath'] = os.path.join(pars['slurm_files_dir_path'], submission_filename)
            
            activate_slurm_submit(pars)
        
        print('\nThe PE config files in {configs_path} are running in {n_jobs} detached slurm jobs, within {conda_env} conda environment. Good luck!\n'.format(configs_path = pars['config_files_dir_path'], conda_env=pars['conda_env'], n_jobs=len(config_filenames)))

# ---------------------------------------------------------------------- #

pars = {
    'conda_env'    : 'icarogw',
    'user_mail'    : '',
    'slurm_nodes'  : 1,
    'slurm_cpus'   : 8,
    'slurm_memory' : 8,
    'slurm_time'   : {'days': 0, 'hours': 8, 'minutes': 0},
}
pars['slurm_python_path']     = '/sps/virgo/USERS/tbertheas/conda/envs/{conda_env}/bin/python'.format(conda_env=pars['conda_env'])
pars['slurm_executable_file'] = '/sps/virgo/USERS/tbertheas/icarogw_pipeline/icarogw_pipeline/parameter_estimation/bilby_pipeline.py'

# MAIN INPUT: Set the directory where the population is stored
pars['population_dir_path']   = 'path_to_population_directory'
# OR
# MAIN INPUT: Set the directory where the config files are stored
pars['config_files_dir_path'] = 'path_to_config_files_directory'

# ---------------------------------------------------------------------- #

if __name__ == '__main__':
    main()
