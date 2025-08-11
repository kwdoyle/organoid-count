#!/bin/bash

toprocess=$1
maindir=$2
basedir='./out/'$toprocess

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
    python ./python_scripts/run_cell_count_set_sections.py "$fl" "$savedir" "$toprocess"
done
