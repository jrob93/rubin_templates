#!/bin/bash

# pass memory per task from the command line
# check if arg is passed via command line
if [ $# -eq 0 ]
  then
    mem="2GB"
else
  mem=$1
fi
echo $mem

# pass a db file keyword from the command line
dbname=$2
# check if db name is passed via command line
if [ $# -eq 0 ]
  then
    dbname=""
else
  dbname=$2
fi
echo $dbname

## declare an array variable
declare -a pop_arr=("mba_5k" "l7_5k" "granvik_5k" "granvik_pha_5k" "occ_rmax5_5k" "occ_rmax20_5k")

nodes=1
ntasks=1
out="test.slurm"

# ensure output file is empty
touch $out
echo -n > $out
 
# write some filler SBATCH stuff
echo "$(cat slurm_header.txt)" >> $out
echo -e "cd $PWD\n" >> $out

# loop over all db files in this directory
for pop in "${pop_arr[@]}"; do
  for f in ss_script*$dbname*$pop.sh; do
    srun="srun --export=ALL --ntasks=$ntasks --nodes=$nodes --mem-per-cpu=$mem ./$f &"
    #srun="srun --export=ALL --ntasks=$ntasks --nodes=$nodes --mem-per-cpu=$mem echo \"hello world\" &"
    echo $srun >> $out
    #echo $srun
    echo "echo \"launch $f\"" >> $out
    wait
  done
done

echo "wait" >> $out

cat $out
