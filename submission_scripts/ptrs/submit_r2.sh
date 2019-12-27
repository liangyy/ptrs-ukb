if [[ ! -d logs ]]
then
  mkdir logs
fi

for i in `seq 1 17`
do
  qsub -v NUM=$i -N R2-$i run_r2.qsub
done