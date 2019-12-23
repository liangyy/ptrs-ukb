if [[ ! -d logs ]]
then
  mkdir -p logs
fi

for i in `seq 1 22`
do 
  qsub -v CHRNUM=$i -N cleanup-$i run_cleanup.qsub
done
