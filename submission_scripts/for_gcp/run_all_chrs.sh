# sublist name: up, middle, bottom
SUB=$1
hailctl dataproc start all-chrs-$SUB \
  --master-boot-disk-size 100 \
  -p 60 \
  --worker-machine-type n1-highmem-16 \
  --master-machine-type n1-highmem-16 \
  --preemptible-worker-boot-disk-size 16 \
  > all-chrs-$SUB.log 2>&1 

hailctl dataproc submit all-chrs-$SUB \
  --pyfiles ../../code/my_hail_helper.py,../../code/gwas_helper.py \
  ../../code/gwas_on_subset_ht.py \
    --bgen-path gs://ukb_data/genotypes/v3/ukb_imp_chr{chr_num}_v3.bgen \
    --bgen-sample gs://ukb_data/genotypes/ukb19526_imp_chr1_v3_s487395.sample \
    --bgen-index gs://ukb_data/bgen_idx/ukb_imp_chr{chr_num}_v3.bgen.idx2 \
    --hail-block-size 128 \
    --output-filename gs://ptrs-ukb/results/gwas_on_subsets/all_chrs_$SUB \
    --variant-ht gs://ukb_data/variant_qc/imp_all.ht \
    --variant-qc-yaml gs://ptrs-ukb/results/misc/gwas_full_variant_qc.yaml \
    --pheno-covar-path gs://ptrs-ukb/results/pheno_and_covar/all_subsets.{subset_name}.ht \
    --subset-yaml gs://ptrs-ukb/results/misc/gwas_full_subset_list_$SUB.yaml \
    --google-cloud-project ukb-im \
    --hail-log /home/yanyu018/hail-all_chrs_$SUB.log >> all-chrs-$SUB.log 2>&1 

# hailctl dataproc stop all-chrs-$SUB

screen -X kill