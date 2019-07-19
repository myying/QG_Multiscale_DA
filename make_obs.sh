#!/bin/bash
export CONFIG=/glade/work/mying/qgmodel_enkf/config/$1/noda
. $CONFIG

mkdir -p $workdir/obs
./generate_obs.py $workdir $kmax $nz $num_cycle $obs_thin
