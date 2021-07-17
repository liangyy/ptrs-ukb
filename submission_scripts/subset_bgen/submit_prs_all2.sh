if [[ ! -f sample_list.prs_all2.txt ]]
then
  cat /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ptrs-ukb/output/data_split/British-test-1.txt | tail -n +2 > sample_list.prs_all2.txt
  cat /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ptrs-ukb/output/data_split/British-validation-1.txt | tail -n +2 >> sample_list.prs_all2.txt
  cat /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ptrs-ukb/output/data_split/African.txt | tail -n +2 >> sample_list.prs_all2.txt
  cat /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ptrs-ukb/output/data_split/Chinese.txt | tail -n +2 >> sample_list.prs_all2.txt
  cat /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ptrs-ukb/output/new_data_split/Indian.txt | tail -n +2 >> sample_list.prs_all2.txt
  cat /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ptrs-ukb/output/new2_data_split/Caribbean.txt | tail -n +2 >> sample_list.prs_all2.txt
fi

for i in `seq 1 22`
do
  echo qsub -v NUM=$i run_subset_bgen_prs_all2.qsub
done
