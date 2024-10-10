import os, sys

template = """executable = {executable}
universe = vanilla
arguments = "{script} --config-file {config}"
request_cpus = {cpus}
request_memory = {memory} GB
request_disk = {disk} GB
log = {condor_path}/{log_name}.log
error = {condor_path}/{log_name}.err
accounting_group = ligo.prod.o4.cbc.hubble.icarogw
accounting_group_user = {user_name}
getenv = True
queue"""


def activate_condor_submit(config_name):

    subfile = '{}/submit_icarogw_{}.sub'.format(sub_path, config_name.split('/')[-1].split('.ini')[0].split('config_')[-1])
    sys.stderr.write('generating {}\n'.format(subfile))

    with open(subfile,'w') as f:

        submission_command = template.format(executable  = condor_executable_path,
                                             script      = condor_executable_file,
                                             config      = config_name,
                                             condor_path = con_path,
                                             log_name    = config_name.split('/')[-1].split('.ini')[0].split('config_')[-1],
                                             cpus        = condor_threads,
                                             memory      = condor_memory,
                                             disk        = condor_disk,
                                             user_name   = user_name)
                                             
        f.write(submission_command)
    sys.stderr.write('submitting {}\n'.format(subfile))
    os.system('condor_submit -batch-name {} {}'.format(subfile.split('/')[-1].split('submit_')[-1].split('.sub')[0], subfile))

# ---------------------------------------------------------------------- #
user_name      = 'vasco.gennari'
condor_memory  = 5
condor_disk    = 1
condor_threads = 10
condor_executable_path = '/home/vasco.gennari/.conda/envs/icarogw_env/bin/python'
condor_executable_file = '/home/vasco.gennari/icarogw_pipeline/icarogw_pipeline/icarogw_runner.py'

# Set the specific directory for the runs
directory    = '/home/vasco.gennari/icarogw_pipeline/config_files'
subdirectory = 'tmp'
# ---------------------------------------------------------------------- #

sub_path = os.path.join(directory, 'submission_files')
con_path = os.path.join(directory, 'condor_files')
if not os.path.exists(sub_path): os.makedirs(sub_path)
if not os.path.exists(con_path): os.makedirs(con_path)

if not (subdirectory == ''): final_path = os.path.join(directory, subdirectory)
else:                        final_path = directory
configs_path = os.path.join(os.getcwd(), final_path)
config_list  = os.listdir(configs_path)

print('')
for config in config_list:
    config_path = os.path.join(configs_path, config)
    activate_condor_submit(config_path)

print('\nThe config files in {configs_path} are running in detached condor jobs. Good luck!\n'.format(configs_path = configs_path))
