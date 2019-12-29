if [[ ! -d logs ]]
then
  mkdir -p logs
fi

for i in `seq 1 17`
do 
  qsub -v SUBSET=subset$i -N merge-$i run_merge.qsub
done

