import os, sys

template = """executable = {executable}
universe = vanilla
arguments = "--config-file {config}"
request_cpus = {cpus}
request_memory = {memory} GB
request_disk = {disk} GB
log = condor_files/{log_name}.log
error = condor_files/{log_name}.err
accounting_group = ligo.dev.o4.cbc.testgr.tiger
getenv = True
queue"""


def activate_condor_submit(config_name):

    subfile = 'submission_files/submit_pyring_{}.sub'.format(config_name.split('/')[-1].split('.ini')[0].split('config_')[-1])
    sys.stderr.write('generating {}\n'.format(subfile))

    with open(subfile,'w') as f:

        submission_command = template.format(executable = os.path.join(condor_executable_path,'pyRing'),
                                             config     = config_name,
                                             log_name   = config_name.split('/')[-1].split('.ini')[0].split('config_')[-1],
                                             cpus       = condor_threads,
                                             memory     = condor_memory,
                                             disk       = condor_disk)
                                             
        f.write(submission_command)
    sys.stderr.write('submitting {}\n'.format(subfile))
    os.system('condor_submit -batch-name {} {}'.format(subfile.split('/')[-1].split('submit_pyring_')[-1].split('.sub')[0], subfile))

# ---------------------------------------------------------------------- #
condor_memory  = 32
condor_threads = 1
condor_disk    = 1
condor_executable_path = '/home/vasco.gennari/.conda/envs/pyring_env/bin'

# Set the specific directory for the runs
directory    = 'SXS_injections'
subdirectory = 'tmp'
# ---------------------------------------------------------------------- #

if not os.path.exists('submission_files'): os.makedirs('submission_files')
if not os.path.exists('condor_files'):     os.makedirs('condor_files')

if not (subdirectory == ''): final_path = os.path.join(directory, subdirectory)
else:                        final_path = directory
configs_path = os.path.join(os.getcwd(), final_path)
config_list  = os.listdir(configs_path)

for config in config_list:
    config_path = os.path.join(configs_path, config)
    activate_condor_submit(config_path)

print('\nThe config files in {configs_path} are running in detached condor jobs. Good luck!'.format(configs_path=configs_path))
