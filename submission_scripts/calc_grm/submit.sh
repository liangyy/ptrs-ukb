for i in `seq 1 22`
do
  qsub -v NUM=$i run.qsub
done