import os

def generate_htcondor_job_scripts():
    tags = [None, 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
    #tags = ['J']
    opt_param_list =  ['kappa_OC', 'kappa_IM', 'kappa_OM']
    #opt_param_list =  ['kappa_IM']
    scaling_flag = [True, False]

    # tags = [None]
    # opt_param_list =  ['kappa_OC']
    # scaling_flag = [True]

    # tags = [None]
    # opt_param_list = ['t13', 'mAtm', 't23', 'dcp']
    # scaling_flag = [False]

    if not os.path.exists('job_scripts'):
        os.mkdir('job_scripts')

    if not os.path.exists('HTCondor_OutErrLog_files'):
        os.mkdir('HTCondor_OutErrLog_files')

    job_index = 0
    for p in opt_param_list:
        for tag in tags:
            for flag in scaling_flag:
                if flag or ((flag == False) and (tag is None)):
                    job_index += 1
                    sub_script_content = f"""notification            = Never
    universe                = vanilla
    executable              = /afs/cern.ch/work/c/cjesus/private/EarthTomo/lxplus_jobs/job_scripts/job_{job_index}.sh
    output                  = job_{job_index}.$(ClusterId).$(ProcId).out
    error                   = job_{job_index}.$(ClusterId).$(ProcId).error
    log                     = job_{job_index}.$(ClusterId).$(ProcId).log
    getenv                  = True
    should_transfer_files   = NO
    initialdir              = /afs/cern.ch/work/c/cjesus/private/EarthTomo/lxplus_jobs/HTCondor_OutErrLog_files
    priority                = 20
    request_cpus            = 2
    request_memory          = 8 GB
    request_disk            = 8 GB
    +JobFlavour             = "tomorrow"
    +MaxRuntime             = 57600
    queue 1
    """
                    with open(f'job_scripts/job_{job_index}.sub', 'w') as f:
                        f.write(sub_script_content)

                    sh_script_content = f"""#!/bin/bash

conda activate /afs/cern.ch/work/c/cjesus/miniforge3/envs/EarthTomo

python - <<END
import time
from app import optManager

opt = optManager(scale_to_HK={flag}, tag={repr(tag)})
stime = time.time()
outname = opt.calc_1D_profile_DeltaChi2('{p}')
print(outname)
print(f'Execution time: {{time.time()-stime:.0f}} seconds.')
END
"""
                    with open(f'job_scripts/job_{job_index}.sh', 'w') as f:
                        f.write(sh_script_content)

    print(f"Generated {job_index} HTCondor job scripts in the 'job_scripts' directory.")

if __name__ == '__main__':
    generate_htcondor_job_scripts()