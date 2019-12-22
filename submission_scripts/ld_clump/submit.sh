for i in `seq 1 22`
do 
  qsub -v CHRNUM=$i -N cleanup-$i run.qsub
done
