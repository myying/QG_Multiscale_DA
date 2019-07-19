#!/bin/bash
source ~/.bashrc

export CONFIG=/glade/work/mying/qgmodel_enkf/config/$1
. $CONFIG

PBS_NP=`cat $PBS_NODEFILE |wc -l`
if [ $casename == "sl" ]; then
  export casename="sl$localize_cutoff"
  if [ $localize_cutoff -eq 256 ]; then
    export localize=0
  fi
fi

mkdir -p $workdir/$casename
cd $workdir/$casename

#first guess is truth+perturbation
#$homedir/add_perturb.sh initial_condition first_guess 3 10 0.8

#get initial prior ensemble from noda case
if [ ! -f current_cycle ]; then
	for m in `seq 1 $nens`; do
		mid=`printf %4.4i $m`
		mkdir -p $mid
		cp $workdir/noda/$mid/f_`printf %5.5i $((spinup_cycle+1))`.bin $mid/.
	done
	echo $((spinup_cycle+1)) > current_cycle
fi

current_cycle=`cat current_cycle`

for n in `seq $current_cycle $num_cycle`; do
echo $n
  if $run_filter && [ $n -gt $spinup_cycle ] && [ $(((n+$obs_interval-1)%$obs_interval)) -eq 0 ]; then
    cd $workdir/$casename

    nt=64
    t=0
    for m in `seq 1 $nens`; do
      rm -f `printf %4.4i $m`/`printf %5.5i $n`.bin &
      t=$((t+1))
      if [ $t -gt $nt ]; then t=0; wait; fi
    done
    wait

    ####Run filter
    $homedir/filter.py $workdir $casename $kmax $nz $nens $n $run_displacement 0 >> filter.log 2>&1

		for m in `seq 1 $nens`; do
			mem=`printf %4.4i $m`
			if [ ! -f $mem/`printf %5.5i $n`.bin ]; then
        echo "filter failed at cycle $n $mem/`printf %5.5i $n`.bin"; exit;
      fi
		done
  else
    nt=$PBS_NP; t=0;
		for m in `seq 1 $nens`; do
			mem=`printf %4.4i $m`
      cd $mem
      ln -fs f_`printf %5.5i $n`.bin `printf %5.5i $n`.bin &
      t=$((t+1))
      if [ $t -gt $nt ]; then t=0; wait; fi
      cd ..
    done
  fi

  nt=$PBS_NP; t=0;
  for m in `seq 1 $nens`; do
    mem=`printf %4.4i $m`
    cd $mem
    rm -f output.bin
    cp -L `printf %5.5i $n`.bin input.bin
    rm -f restart.nml
    $homedir/namelist_input.sh $n 1 > input.nml #same idum as truth (same random forcing)
    $codedir/$qgexe . >& /dev/null &
    t=$((t+1))
    if [ $t -gt $nt ]; then t=0; wait; fi
    cd ..
  done
	wait

  for m in `seq 1 $nens`; do
  	mem=`printf %4.4i $m`
    mv $mem/output.bin $mem/f_`printf %5.5i $((n+1))`.bin
		if [ ! -f $mem/f_`printf %5.5i $((n+1))`.bin ]; then
      echo "filter failed at cycle $n $mem/`printf %5.5i $((n+1))`.bin"; exit;
    fi
	done

	echo $((n+1)) > current_cycle
done
