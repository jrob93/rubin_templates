#!/bin/bash

# pass population from the command line
pop=$1
# check if file name is passed via command line
if [ $# -eq 0 ]
  then
    pop="mba_5k"
else
  pop=$1
fi
echo $pop

# create a tmp file to store the header inforamtion for cuillin (conda etc)
echo $PWD
cp ss_script_header.txt ss_script_header.tmp
echo "cd $PWD" >> ss_script_header.tmp 

# loop over all db files in this directory
for f in *.db; do
  echo "db -> $f, pop -> $pop"
  
  # get the new name to save the script 
  arrIN=(${f//.db/ })
  f2=ss_script_${arrIN[0]}_$pop.sh
  echo "out -> $f2"

  # command to generate ss_script.sh
  cmd1="generate_ss --db $f --pop $pop"
  echo "run -> $cmd1"
  $cmd1

  # command to rename ss_script.sh
  cmd2="mv ss_script.sh $f2"
  echo "run -> $cmd2"
  $cmd2

  # combine the header and ss_script file
  # use temporary file out.txt
  cat ss_script_header.tmp $f2 >> out.txt
  mv out.txt $f2

  # make the script executable
  chmod +x $f2

  wait
done

# remove the tmp file
rm ss_script_header.tmp
