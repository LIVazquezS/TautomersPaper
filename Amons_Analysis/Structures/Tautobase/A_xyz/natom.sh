#!/bin/bash

for i in {1..1259}; do

nc=$(sed -n '/C/=' a_${i}.xyz | wc -l)
no=$(sed -n '/O/=' a_${i}.xyz | wc -l)
nn=$(sed -n '/N/=' a_${i}.xyz | wc -l)
nHeavy=$(echo "$nc + $nn + $no" | bc -l)

echo $i $nc >> nc_tauto.txt
echo $i $nn >> nn_tauto.txt
echo $i $no >> no_tauto.txt
echo $i $nHeavy >> nheavy_tauto.txt

done 
