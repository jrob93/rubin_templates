# rubin_templates

run template_metrics_slurm.py to generate the slurm scripts which calls template_metrics.py

template_metrics.py has several options to run the template generation metrics

Run analyse_no_template_visits on cuillin to populate the analysis_metric_bundles dir
rsync -nav -e ssh cuillin:/home/jrobinson/rubin_templates/analysis_metric_bundles* .

Run the ss_metrics analysis

Download all ss_metric npz files to /Volumes/Nargothrond/rubin_templates/ss_metrics

Download all resultsDB_sqlite.db files to /Users/jrobinson/rubin_templates/ss_metrics
rsync -nav -e ssh cuillin:/home/jrobinson/rubin_templates/ss_metrics/*baseline_v3_3*noTwi* ss_metrics --exclude="*.out" --exclude="*baseline*.db" --exclude="*.slurm" --exclude="*.sh" --exclude="*.png" --exclude="*.npz" --exclude="*.pdf" --exclude="*.txt"
rsync -nav -e ssh cuillin:/home/jrobinson/rubin_templates/ss_metrics/first_year*noDD_noTwi*frac90_ss ss_metrics --exclude="*.out" --exclude="*baseline*.db" --exclude="*.slurm" --exclude="*.sh" --exclude="*.png" --exclude="*.npz" 
--exclude="*.pdf" --exclude="*.txt"
