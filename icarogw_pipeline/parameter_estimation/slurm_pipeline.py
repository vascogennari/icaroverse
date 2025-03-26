import os, sys

#SBATCH --tasks-per-node=1

template = """#!/bin/sh

#SBATCH --job-name={name}
#SBATCH --output={event_dir_path}/slurm_logs/output_{name}.out
#SBATCH --error={event_dir_path}/slurm_logs/error_{name}.err
#SBATCH --nodes={nodes}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory}G
#SBATCH --partition=htc
#SBATCH --time={time}
#SBATCH --account=virgo
#SBATCH --licenses=sps
#SBATCH --mail-user={user_mail}
#SBATCH --mail-type=ALL

module load conda
conda activate {conda_env}
{executable} {script} --config-file {config}
"""

def activate_slurm_submit(event_dir_path, population_dir_path):

    submission_filename = os.path.join(event_dir_path, "slurm_submission_file.sh")
    sys.stderr.write('generating {}\n'.format(submission_filename))

    config_filename = os.path.join(event_dir_path, "config_for_bilby_PE.ini")

    with open(submission_filename, 'w') as f:
        submission_command = template.format(
            name            = '_'.join([os.path.basename(population_dir_path), os.path.basename(event_dir_path)]),
            event_dir_path  = event_dir_path,
            nodes           = slurm_nodes,
            cpus            = slurm_cpus,
            memory          = slurm_memory,
            time            = '{}-{}:{}:00'.format(slurm_time['days'], slurm_time['hours'], slurm_time['minutes']),
            user_mail       = user_mail,
            conda_env       = conda_env,
            executable      = slurm_executable_path,
            script          = slurm_executable_file,
            config          = config_filename)
                                             
        f.write(submission_command)
    
    sys.stderr.write('submitting {}\n\n'.format(submission_filename))
    os.system('sbatch {}'.format(submission_filename))

# ---------------------------------------------------------------------- #
conda_env    = 'in2_env'
user_mail    = 'tom.bertheas@l2it.in2p3.fr'
slurm_nodes  = 1
slurm_cpus   = 10
slurm_memory = 5
slurm_time   = {'days': 2, 'hours': 0, 'minutes': 0}
slurm_executable_path = '/pbs/home/t/tbertheas/.conda/envs/{conda_env}/bin/python'.format(conda_env=conda_env)
slurm_executable_file = '/sps/virgo/USERS/tbertheas/icarogw_pipeline/icarogw_pipeline/parameter_estimation/bilby_pipeline.py'

# MAIN INPUT: Set the directory where the population is stored
population_dir_path    = '<population_dir_path>'
# ---------------------------------------------------------------------- #

population_PE_dir_path =  os.path.join(population_dir_path, "parameter_estimation")
if not os.path.exists(population_PE_dir_path):
    raise FileNotFoundError("No `parameter_estimation` directory found in the population directory. Please make sure the `population_directory` path is correct and/or the corresponding population has been processed with `generate_configs.py`")

events_dir_list = os.listdir(population_PE_dir_path)

print('')
for event_dir_name in events_dir_list:
    event_dir_path = os.path.join(population_PE_dir_path, event_dir_name)
    activate_slurm_submit(event_dir_path, population_dir_path)

print('\nThe PE config files in {configs_path} are running in detached slurm jobs, within {conda_env} conda environment. Good luck!\n'.format(configs_path = population_PE_dir_path, conda_env=conda_env))
