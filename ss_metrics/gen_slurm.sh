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

# also pass the population
if [ $# -eq 0 ]
  then
    pop="mba_5k"
else
  pop=$2
fi
echo $pop

# pass a db file keyword from the command line
dbname=$3
# check if db name is passed via command line
if [ $# -eq 0 ]
  then
    dbname=""
else
  dbname=$3
fi
echo $dbname

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
for f in ss_script*$dbname*$pop.sh; do
  srun="srun --export=ALL --ntasks=$ntasks --nodes=$nodes --mem-per-cpu=$mem ./$f &"
  #srun="srun --export=ALL --ntasks=$ntasks --nodes=$nodes --mem-per-cpu=$mem echo \"hello world\" &"
  echo $srun >> $out
  #echo $srun
  echo "echo \"launch $f\"" >> $out
  wait
done

echo "wait" >> $out

cat $out
