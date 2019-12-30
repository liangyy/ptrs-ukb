NAME=prs-makeup
JOBNAME=$NAME

hailctl dataproc start $JOBNAME \
  --master-boot-disk-size 100 \
  -p 40 \
  --zone us-central1-b \
  --worker-machine-type n1-highmem-16 \
  --master-machine-type n1-highmem-16 \
  --preemptible-worker-boot-disk-size 16 \
  > $JOBNAME.log 2>&1 

# # hailctl dataproc start $JOBNAME --master-boot-disk-size 60 -p 20 --preemptible-worker-boot-disk-size 16 > $JOBNAME.log 2>&1 

PREFIX=3
i=10_makeup
hailctl dataproc submit $JOBNAME \
  --pyfiles ../../code/prs_helper.py,../../code/gwas_helper.py \
  ../../code/ld_clump_prs.py \
    --mode skip_subset \
    --dont-overwrite \
    --bgen-path gs://ukb_data/genotypes/v3/ukb_imp_chr{chr_num}_v3.bgen \
    --bgen-sample gs://ukb_data/genotypes/ukb19526_imp_chr1_v3_s487395.sample \
    --bgen-index gs://ukb_data/bgen_idx/ukb_imp_chr{chr_num}_v3.bgen.idx2 \
    --pval-thresholds '5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,0.05,0.1,0.5,1' \
    --variant-pool gs://ptrs-ukb/results/gwas_on_subsets/in_tsv_format/gwas_in_tsv_subset3_x_height.tsv \
    --google-cloud-project ukb-im \
    --subset-and-gwas gs://ptrs-ukb/results/misc/prs_gcp_$i.yaml \
    --mt-prefix gs://ptrs-ukb/results/prs/presubset_$PREFIX \
    --output-prefix gs://ptrs-ukb/prs/prs_gcp \
    --hail-log gs://ptrs-ukb/tmp-output/hail_prs_gcp_$i.log \
    >> $JOBNAME-$i.log 2>&1 


PREFIX=4
i=14
hailctl dataproc submit $JOBNAME \
  --pyfiles ../../code/prs_helper.py,../../code/gwas_helper.py \
  ../../code/ld_clump_prs.py \
    --mode skip_subset \
    --dont-overwrite \
    --bgen-path gs://ukb_data/genotypes/v3/ukb_imp_chr{chr_num}_v3.bgen \
    --bgen-sample gs://ukb_data/genotypes/ukb19526_imp_chr1_v3_s487395.sample \
    --bgen-index gs://ukb_data/bgen_idx/ukb_imp_chr{chr_num}_v3.bgen.idx2 \
    --pval-thresholds '5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,0.05,0.1,0.5,1' \
    --variant-pool gs://ptrs-ukb/results/gwas_on_subsets/in_tsv_format/gwas_in_tsv_subset3_x_height.tsv \
    --google-cloud-project ukb-im \
    --subset-and-gwas gs://ptrs-ukb/results/misc/prs_gcp_$i.yaml \
    --mt-prefix gs://ptrs-ukb/results/prs/presubset_$PREFIX \
    --output-prefix gs://ptrs-ukb/prs/prs_gcp \
    --hail-log gs://ptrs-ukb/tmp-output/hail_prs_gcp_$i.log \
    >> $JOBNAME-$i.log 2>&1 

  
# --gwas-ht gs://ptrs-ukb/results/gwas_on_subsets/all_chrs_up.ht \
hailctl dataproc stop $JOBNAME

