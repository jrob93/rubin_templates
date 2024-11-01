import numpy as np

# base_cmd = "python -u template_metrics.py"
# jobname_suffix = ""

# base_cmd = "python -u template_metrics_n_visits_4.py"
# jobname_suffix = "_n_visits_4"

base_cmd = "python template_metrics_n_visits_4_override-g.py"
jobname_suffix = "_n_visits_4_override_g"

"""
# db_list = ["baseline_v3.0_10yrs.db","ender_a1_v3.1_10yrs.db","baseline_v3.2_10yrs.db"]
db_list = ["ender_a1_v3.1_10yrs.db"]
nside_list = [256]
tscales = [3,7,14,28]
# runs = ["baseline", "metrics", "pairs", "visits"]
runs = ["metrics", "pairs", "visits"]

db_list = ["baseline_v3.3_10yrs.db"]
nside_list = [256]
tscales = [7,14]
runs = ["metrics", "visits"]
"""

db_list = ["baseline_v3.3_10yrs.db"]
nside_list = [256]
tscales = [3,28]
runs = ["metrics","visits"]

# slurm setup
mem_per_task = 2
srun_cmd = "srun --export=ALL --ntasks=1 --nodes=1 --mem-per-cpu={}GB --exclusive".format(mem_per_task)
max_hours = 5*24
run_dir = "/home/jrobinson/rubin_templates"

for t in tscales:
    for d in db_list:
        for n in nside_list:
            for run in runs:

                count = 0
                cmd_list = []
                job_name = "{}_{}_{}_{}{}".format("_".join(d.split(".")),n,t,run,jobname_suffix)

                # out_file = "{}_{}_{}_{}".format("_".join(d.split(".")),n,t,run)
                out_file = job_name
                cmd = "{} -d {} -n {} -t {} --{} > {}.out &\necho \"launch {}\"".format(base_cmd,d,n,t,run,out_file,out_file)
                cmd = srun_cmd + " " + cmd
                cmd_list.append(cmd)
                count+=1

                print("\n{} runs\n".format(count))

                # make the slurm script
                slurm_name = "{}.slurm".format(job_name)

                if t==3:
                    _max_hours = max_hours*2
                else:
                    _max_hours = max_hours

                slurm_script="""#!/bin/bash \

#SBATCH --job-name={}
#SBATCH --nodes=1
#SBATCH --ntasks={}
#SBATCH --cpus-per-task=1
#SBATCH --time={}:00:00
#SBATCH --mem={}GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.robinson@ed.ac.uk

# . /usr/local/anaconda/3.9/etc/profile.d/conda.sh
# conda activate rubin
source /home/jrobinson/miniconda3/etc/profile.d/conda.sh
conda activate rubin2
cd {}

{}\nwait""".format(job_name,count,_max_hours,count*mem_per_task,run_dir,"\n".join(cmd_list))

                print(slurm_script)

                with open(slurm_name,"w") as f:
                    f.write(slurm_script)
