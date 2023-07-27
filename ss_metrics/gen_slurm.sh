#!/bin/bash

# pass memory per task from the command line
# check if arg is passed via command line
if [ $# -eq 0 ]
  then
    mem="2GB"
else
  mem=$1
fi

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
for f in ss_script*.sh; do
  srun="srun --export=ALL --ntasks=$ntasks --nodes=$nodes --mem-per-cpu=$mem ./$f &"
  #srun="srun --export=ALL --ntasks=$ntasks --nodes=$nodes --mem-per-cpu=$mem echo \"hello world\" &"
  echo $srun >> $out
  #echo $srun
  echo "echo \"launch $f\"" >> $out
  wait
done

echo "wait" >> $out

cat $out
