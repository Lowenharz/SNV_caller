#~/bin/bash

file_dir=$1 

if [ ! -d on_target_craft ]; then 
    mkdir -p on_target_craft
fi

for file in $1/fold.change.*.txt; do 
    if [ ! -f $file ]; then 
        echo $file not exist. Make sure you are entering the correct directory. 
        exit
    fi
    base_name=$(basename $file)
    awk '$5>20' $file |cat >> on_target_craft/on_target_$base_name
done 