#!/bin/bash
. /usr/local/anaconda/3.9/etc/profile.d/conda.sh
conda activate rubin
cd /home/jrobinson/rubin_templates/ss_metrics/
make_lsst_obs --simulation_db first_year_baseline_v3_0_10yrs.db --orbit_file /home/jrobinson/rubin_sim_data/orbits/mba_5k.txt ; run_moving_calc --obs_file first_year_baseline_v3_0_10yrs__mba_5k_obs.txt --simulation_db first_year_baseline_v3_0_10yrs.db --orbit_file /home/jrobinson/rubin_sim_data/orbits/mba_5k.txt --out_dir first_year_baseline_v3_0_10yrs_ss --objtype MBA --start_time 60218.0 ; run_moving_fractions --work_dir first_year_baseline_v3_0_10yrs_ss --metadata MBA --start_time 60218.0
