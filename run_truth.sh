#!/bin/bash
#PBS -A UPSU0023
#PBS -N run_truth
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=16
#PBS -q regular
#PBS -j oe
#PBS -o log
source ~/.bashrc

export CONFIG=/glade/work/mying/qgmodel_enkf/config/$1/noda
. $CONFIG

mkdir -p $workdir/truth
cd $workdir/truth

if [ ! -f current_cycle ]; then
  echo 1 > current_cycle
  cp $workdir/initial_condition.bin input.bin
fi

current_cycle=`cat current_cycle`
for n in `seq $current_cycle $num_cycle`; do
  echo $n

  cp input.bin `printf %5.5i $n`.bin
  rm -f output.bin

  #each cycle step uses a different idum
	idum=$n  #truth model stochastic forcing
	#idum=`echo "(14423*$n+5324)%13431" |bc` #rand error in stoch. forcing
  $homedir/namelist_input.sh $idum > input.nml 
  rm -f restart.nml

  $codedir/$qgexe . >& /dev/null

  if [ ! -f output.bin ]; then 
    echo 'abort'
    exit
  fi
  mv output.bin input.bin
  
  echo $((n+1)) > current_cycle
done
