if [[ ! -d logs ]]
then
  mkdir -p logs
fi

for i in `cat ../../misc/trait_list.txt`
do 
  qsub -v TRAIT=$i -N subset-gwas--$i run_subset_gwas.qsub
done

