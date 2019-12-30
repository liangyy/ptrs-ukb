# ARGS1: OUTDIR

if [[ ! -d logs ]]
then
  mkdir logs
fi

for i in `seq 1 3`
do
  qsub -v NUM=$i,CONFIG=no_name,OUTDIR=$1 -N R2-$i run.qsub
done

for i in `seq 4 17`
do
  qsub -v NUM=$i,CONFIG=with_name,OUTDIR=$1 -N R2-$i run.qsub
done
