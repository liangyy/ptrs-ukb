for i in `seq 1 22`
do
  qsub -v NUM=$i -N CHR$i run.qsub
done
