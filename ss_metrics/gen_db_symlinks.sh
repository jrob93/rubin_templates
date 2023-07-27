#!/bin/bash

# read the file with all dbs for processing
# ls ../*/*.db > visit_db_list.txt

# check if file name is passed via command line 
if [ $# -eq 0 ]
  then
    db_file=visit_db_list.txt
else
  db_file=$1
fi

echo $db_file

while read f1; do
  echo "$f1"

  # get the db name to create the symlink
  IFS='/'; arrf1=($f1); unset IFS;
  f2=${arrf1[${#arrf1[@]}-1]}

  # create and run the symlink command
  cmd="ln -s $f1 $f2"
  echo $cmd
  $cmd

done <$db_file
