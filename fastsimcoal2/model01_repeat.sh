#!/bin/sh


PREFIX="model01"

for ITER in {1..100}
do

MY_WORKDIR=run${ITER}

   mkdir $MY_WORKDIR
   cp ${PREFIX}.tpl ${PREFIX}.est ${PREFIX}_MSFS.obs $MY_WORKDIR
   cd $MY_WORKDIR
   /home/user/Software/fsc26_linux64/fsc26 -t ${PREFIX}.tpl -e ${PREFIX}.est  -n 1000000 -C 10 -M -L 40 -c 20   -m --multiSFS -q 
   cd ..
done
