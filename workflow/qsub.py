#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

# do something useful with the threads
threads = job_properties['threads']
resources = int(job_properties['resources']['mem_gb'])

if resources < 50:
    os.system("qsub -cwd -l vf={resources},num_proc={threads} -P P21Z25400N0107 -q st.q {script}".format(resources = resources, threads=threads, script=jobscript))
else:
    os.system("qsub -cwd -l vf={resources},num_proc={threads} -P st_supermem -q st_supermem.q {script}".format(resources = resources, threads=threads, script=jobscript))
