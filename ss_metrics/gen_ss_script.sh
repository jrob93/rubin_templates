#!/bin/bash

# pass population from the command line
pop=$1

# set up the conda env
. /usr/local/anaconda/3.9/etc/profile.d/conda.sh
conda activate rubin

# loop over all db files
for f in *.db; do
  echo "db -> $f, pop -> $pop"
  
  # get the new name to save the script 
  arrIN=(${f//.db/ })
  f2=ss_script_${arrIN[0]}_$pop.sh
  echo "out -> $f2"

  # command to generate ss_script.sh
  cmd1="generate_ss --db $f --pop $pop"
  echo "run -> $cmd1"
  #$cmd1

  # command to rename ss_script.sh
  cmd2="mv ss_script.sh $f2"
  echo "run -> $cmd2"
  #$cmd2

  wait
done
