import os, sys
from argparse import ArgumentParser, ArgumentTypeError

#SBATCH --tasks-per-node=1

template = """#!/bin/sh

#SBATCH --job-name={name}
#SBATCH --output={slurm_path}/output_logs/output_{name}.out
#SBATCH --error={slurm_path}/error_logs/error_{name}.err
#SBATCH --nodes={nodes}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={memory}G
#SBATCH --partition={partition}
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

def activate_slurm_submit(config_name):

    subfile = '{}/submit_icarogw_{}.sh'.format(sub_path, config_name.split('/')[-1].split('.ini')[0].split('config_')[-1])
    sys.stderr.write('generating {}\n'.format(subfile))

    with open(subfile,'w') as f:
        if user_mail != '':
            email_option = email_option_template.format(user_mail=user_mail)
        else:
            email_option = ''
        submission_command = template.format(
            name         = config_name.split('/')[-1].split('.ini')[0].split('config_')[-1],
            slurm_path   = slurm_path,
            nodes        = slurm_nodes,
            cpus         = slurm_cpus,
            memory       = slurm_memory,
            partition    = slurm_partition,
            time         = '{}-{}:{}:00'.format(slurm_time['days'], slurm_time['hours'], slurm_time['minutes']),
            email_option = email_option,
            conda_env    = conda_env,
            executable   = slurm_executable_path,
            script       = slurm_executable_file,
            config       = config_name
        )

        f.write(submission_command)

    sys.stderr.write('submitting {}\n\n'.format(subfile))
    os.system('sbatch {}'.format(subfile))

parser = ArgumentParser()
parser.add_argument("-g", "--gpu", nargs="?", const='any', default=None, type=str, help="If provided, will request a GPU job, with specified GPU model (available : v100, h100). Otherwise, will request a classical CPU job")
args = parser.parse_args()

# ---------------------------------------------------------------------- #
conda_env    = ''
user_mail    = ''
slurm_nodes  = 1
slurm_cpus   = 4
slurm_gpus   = 1
slurm_memory = 8
slurm_time   = {'days': 7, 'hours': 0, 'minutes': 0}
slurm_executable_path = '/sps/virgo/USERS/vgennari/conda/envs/{conda_env}/bin/python'.format(conda_env=conda_env)
slurm_executable_file = '/sps/virgo/USERS/vgennari/icaroverse/icaroverse/icaroverse_runner.py'

# Handling request of CPU / GPU jobs (CC-IN2P3)
available = ['v100', 'h100']
if args.gpu is None:
    slurm_partition    = "htc"
elif args.gpu == "any":
    slurm_partition = f"gpu_v100,gpu_h100 \n#SBATCH --gpus={slurm_gpus}"
elif args.gpu in available:
    slurm_partition = f"gpu_{args.gpu} \n#SBATCH --gpus={slurm_gpus}"
else:
    raise ArgumentTypeError(f"Invalid GPU type: {args.gpu}. Choose from {available + ["any"]}.")

# Set the specific directory for the runs
directory    = '/sps/virgo/USERS/vgennari/icaroverse/config_files/GWTC-4'
subdirectory = 'splines'
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
