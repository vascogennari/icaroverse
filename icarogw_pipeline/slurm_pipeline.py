import os, sys

#SBATCH --tasks-per-node=1

template = """#!/bin/sh

#SBATCH --job-name={name}
#SBATCH --output={slurm_path}/output_logs/output_{name}.out
#SBATCH --error={slurm_path}/error_logs/error_{name}.err
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
<<<<<<< HEAD
conda activate {conda_env}
=======
conda activate icarogw_env
>>>>>>> main
{executable} {script} --config-file {config}
"""

def activate_slurm_submit(config_name):

    subfile = '{}/submit_icarogw_{}.sh'.format(sub_path, config_name.split('/')[-1].split('.ini')[0].split('config_')[-1])
    sys.stderr.write('generating {}\n'.format(subfile))

    with open(subfile,'w') as f:
        submission_command = template.format(name       = config_name.split('/')[-1].split('.ini')[0].split('config_')[-1],
                                             slurm_path = slurm_path,
                                             nodes      = slurm_nodes,
                                             cpus       = slurm_cpus,
                                             memory     = slurm_memory,
                                             time       = '{}-{}:{}:00'.format(slurm_time['days'], slurm_time['hours'], slurm_time['minutes']),
                                             user_mail  = user_mail,
                                             conda_env  = conda_env,
                                             executable = slurm_executable_path,
                                             script     = slurm_executable_file,
                                             config     = config_name)
                                             
        f.write(submission_command)
    sys.stderr.write('submitting {}\n\n'.format(subfile))
    os.system('sbatch {}'.format(subfile))

# ---------------------------------------------------------------------- #
conda_env    = 'icarogw_env'
user_mail    = 'vasco.gennari@l2it.in2p3.fr'
slurm_nodes  = 1
slurm_cpus   = 10
slurm_memory = 5
slurm_time   = {'days': 2, 'hours': 0, 'minutes': 0}
slurm_executable_path = '/sps/virgo/USERS/vgennari/conda/envs/{conda_env}/bin/python'.format(conda_env=conda_env)
slurm_executable_file = '/sps/virgo/USERS/vgennari/icarogw_pipeline/icarogw_pipeline/icarogw_runner.py'

# Set the specific directory for the runs
directory    = '/sps/virgo/USERS/vgennari/icarogw_pipeline/config_files'
subdirectory = 'simulation_evolving_population_PRODX'
# ---------------------------------------------------------------------- #

sub_path   = os.path.join(directory, 'submission_files')
slurm_path = os.path.join(directory, 'slurm_files')
if not os.path.exists(sub_path):   os.makedirs(sub_path)
if not os.path.exists(slurm_path): os.makedirs(slurm_path)

if not (subdirectory == ''): final_path = os.path.join(directory, subdirectory)
else:                        final_path = directory
configs_path = os.path.join(os.getcwd(), final_path)
config_list  = os.listdir(configs_path)

print('')
for config in config_list:
    config_path = os.path.join(configs_path, config)
    activate_slurm_submit(config_path)

print('\nThe config files in {configs_path} are running in detached slurm jobs, within {conda_env} conda environment. Good luck!\n'.format(configs_path = configs_path, conda_env=conda_env))
