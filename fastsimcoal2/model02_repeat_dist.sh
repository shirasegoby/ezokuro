#!/bin/sh

PREFIX="model02"

for i in {1..100}
 do
   mkdir run$i
   cp ${PREFIX}_maxL.par  ${PREFIX}_MSFS.obs run$i"/"
   cd run$i
   mv ${PREFIX}_MSFS.obs ${PREFIX}_maxL_MSFS.obs
    /home/ashirase/soft/fsc26   -i ${PREFIX}_maxL.par   -n 1000000 -C 10 -L 40 -c 20 -B 20 -m --multiSFS -q -x
   sed -n '2,3p' ${PREFIX}_maxL/${PREFIX}_maxL.lhoods | cut -f 2 >> ../${PREFIX}.lhoods_distribution
   cd .. 
 done
