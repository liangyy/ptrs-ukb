if [[ ! -f sample_list.British-test-1.txt ]]
then
  cat /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ptrs-ukb/output/data_split/British-test-1.txt | tail -n +2 > sample_list.British-test-1.txt
fi

for i in `seq 1 22`
do
  qsub -v NUM=$i run_subset_bgen_eur_test.qsub
done