#!/bin/bash

export CONFIG=/glade/work/mying/qgmodel_enkf/config/$1/noda
. $CONFIG

mkdir -p $workdir/$casename
for m in `seq 1 $nens`; do
  mkdir $workdir/$casename/`printf %4.4i $m`
done

./generate_initial_ensemble.py $workdir $casename $kmax $nz $nens 0.001 1

echo 1 > $workdir/$casename/current_cycle
