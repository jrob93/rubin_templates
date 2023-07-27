#!/bin/bash

#SBATCH --job-name=ss_script_mba_5k
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.robinson@ed.ac.uk

. /usr/local/anaconda/3.9/etc/profile.d/conda.sh
conda activate rubin
cd /home/jrobinson/rubin_templates/ss_metrics/

srun --export=ALL --ntasks=1 --nodes=1 --mem-per-cpu=2GB --exclusive ./ss_script_first_year_baseline_v3_0_10yrs_mba_5k.sh &
echo "launch 1"
srun --export=ALL --ntasks=1 --nodes=1 --mem-per-cpu=2GB --exclusive ./ss_script_first_year_baseline_v3_0_10yrs_visit_cut_t-14d_nside-256_mba_5k.sh &
echo "launch 2"
srun --export=ALL --ntasks=1 --nodes=1 --mem-per-cpu=2GB --exclusive ./ss_script_first_year_baseline_v3_0_10yrs_visit_cut_t-14d_nside-32_mba_5k.sh &
echo "launch 3"
srun --export=ALL --ntasks=1 --nodes=1 --mem-per-cpu=2GB --exclusive ./ss_script_first_year_baseline_v3_0_10yrs_visit_cut_t-28d_nside-256_mba_5k.sh &
echo "launch 4"
srun --export=ALL --ntasks=1 --nodes=1 --mem-per-cpu=2GB --exclusive ./ss_script_first_year_baseline_v3_0_10yrs_visit_cut_t-28d_nside-32_mba_5k.sh &
echo "launch 5"
srun --export=ALL --ntasks=1 --nodes=1 --mem-per-cpu=2GB --exclusive ./ss_script_first_year_baseline_v3_0_10yrs_visit_cut_t-7d_nside-32_mba_5k.sh &
echo "launch 6"
srun --export=ALL --ntasks=1 --nodes=1 --mem-per-cpu=2GB --exclusive ./ss_script_first_year_ender_a1_v3_1_10yrs_mba_5k.sh &
echo "launch 7"
srun --export=ALL --ntasks=1 --nodes=1 --mem-per-cpu=2GB --exclusive ./ss_script_first_year_ender_a1_v3_1_10yrs_visit_cut_t-28d_nside-256_mba_5k.sh &
echo "launch 8"

wait
