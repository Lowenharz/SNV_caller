#!/bin/bash

while read -a line; do
    name=${line[0]}
    col_num=${line[1]}
    # echo $name , $col_num
    awk -v a=$col_num '{print $2,$a}' CCLE_copynumber_byGene_2013-12-03.txt | cat > lite_${name}_CCLE.txt
done < CCLE_cellline_list.txt
