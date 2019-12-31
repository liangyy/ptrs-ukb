# ARGS1: configfile middle name

if [[ ! -d logs ]]
then
  mkdir logs
fi

CONFIG=$1

for i in `seq 1 17`
do
  qsub -v NUM=$i,CONFIG=$CONFIG -N R2-$i run_r2.qsub
done