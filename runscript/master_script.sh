#!/bin/bash

#get the current script directory. DOES not work with simulinks!
#in the case simulinks are used pls specfiy folder by hand
scrip_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#specfiy some paths...
base_dir=$(dirname $scrip_dir)
data_dir=$base_dir/data
results_dir=$base_dir/results

human="Homo sapiens"
mouse="Mus musculus"
chimp="Pan troglodytes"
rat="Rattus norvegicus"
to_exclude="all"

#specify output file paths..
output_2ab=$results_dir/2ab.txt
output_2c=$results_dir/2c.txt
output_3=$results_dir/3.txt

#print some output for logging
echo '#########################################################################'
echo 'Starting meNOGG pipeline'
echo 'Specified base directory: ' $base_dir
echo '#########################################################################'

#switching to the base dir...
cd $base_dir
#check if the data folder is there...
if [ ! -d "$data_dir" ]; then
  echo 'Data directory is not in the expected place.'
  echo 'Are you using simulinks? If so, please specify the path in this script.'
fi
#check if the resuts dir exist, if not make it
if [ ! -d "$results_dir" ]; then
  echo 'Creating results directory'
  mkdir -v $results_dir
fi

#run the first question and print some output...
echo 'Generating results for question 2 a/b)...'
echo 'Output file name: '  $output_2ab
if [ -f "$output_2ab" ]; then
  echo 'File already exists. Skipping...'
else
  python eggnog-db.py -o1 "$human" -o2 "$chimp" -x "$mouse" \
    -i "$data_dir" \
    -o "$output_2ab" -p
fi
echo 'Generating results for question 2 a/b)... - Done'

#run the second question and print some output...
echo 'Generating results for question 2 c)...'
echo 'Output file name: '  $output_2c
if [ -f "$output_2c" ]; then
  echo 'File already exists. Skipping...'
else
  python eggnog-db.py -o1 "$human" -o2 "$chimp" -x "$mouse" \
    -i "$data_dir" \
    -o "$output_2c" -d
fi
echo 'Generating results for question 2 c)... - Done'

#run the third question and print some output...
echo 'Generating results for question 3...'
echo 'Output file name: '  $output_3
if [ -f "$output_3" ]; then
  echo 'File already exists. Skipping...'
else
  python eggnog-db.py -o1 "$mouse" -o2 "$rat" -x "$to_exclude" \
    -i "$data_dir" \
    -o "$output_3" -e
fi
echo 'Generating results for question 3... - Done'
