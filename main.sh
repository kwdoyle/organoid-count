#!/bin/bash

toprocess=$1
maindir=$2
stain_name=$3
param_file=$4

if [ -z "$stain_name" ]; then
    dirnm=$toprocess
else
    dirnm=$stain_name
fi

#basedir='./out/'$toprocess
basedir='./out/'$dirnm

echo Looking for ID folders in: $maindir
echo Will save output in $basedir

for fl in "$maindir"/*; do
    echo $fl
    # the directory name (the sample id)
    id_dir=$(basename "$fl")
    echo $id_dir
    savedir=${basedir}/${id_dir}
    echo Save directory is $savedir

    # run the main analysis script
    python ./python_scripts/run_cell_count_set_sections.py "$fl" "$savedir" "$toprocess" "$stain_name" "$param_file"
done
