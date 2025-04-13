for script in job_scripts/*.sub; do
    condor_submit "$script"
done