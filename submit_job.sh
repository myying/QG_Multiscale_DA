#!/bin/bash
expname=$1
casename=$2
ncpus=$3

nnode=$((($ncpus+63)/64))

cat > tmp.sh << EOF
#!/bin/bash
#PBS -A ???
#PBS -N $casename
#PBS -l walltime=12:00:00
#PBS -l select=$nnode:ncpus=64
#PBS -q regular
#PBS -j oe
#PBS -o log/$expname.$casename
cd /glade/work/mying/qgmodel_enkf
./run_cycle.sh $expname/$casename
EOF

qsub tmp.sh
rm tmp.sh
