if [[ ! -d logs ]]
then
  mkdir -p logs
fi

for i in `cat ../../misc/trait_list.txt`
do 
  qsub -v TRAIT=$i -N ld-clump-$i run_ld_clump.qsub
done
