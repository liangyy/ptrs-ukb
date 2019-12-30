NAME=test
SUBSET=1
JOBNAME=$NAME-subset$SUBSET

hailctl dataproc start $JOBNAME \
  --master-boot-disk-size 100 \
  -p 60 \
  --zone us-central1-b \
  --worker-machine-type n1-highmem-16 \
  --master-machine-type n1-highmem-16 \
  --preemptible-worker-boot-disk-size 16 \
  > $JOBNAME.log 2>&1 

# hailctl dataproc start $JOBNAME --master-boot-disk-size 60 -p 20 --preemptible-worker-boot-disk-size 16 > $JOBNAME.log 2>&1 

hailctl dataproc submit $JOBNAME \
  --pyfiles ../../code/prs_helper.py,../../code/gwas_helper.py \
  ../../code/ld_clump_prs_quick_test.py \
    --pval-thresholds '5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,0.05,0.1,0.5,1' \
    --gwas-tsv gs://ptrs-ukb/results/ld_clump/merged_by_subset.subset1.gwas.tsv \
    --trait-names basophil,bmi,dbp,eosinophil,hb,height,ht,lymphocyte,mch,mchc,mcv,monocyte,neutrophil,platelet,rbc,sbp,wbc \
    --mt-input gs://ptrs-ukb/results/prs/presubset_1_x_subset1.mt \
    --output gs://ptrs-ukb/tmp-output/prs_test_quick.tsv.bgz \
    --hail-log gs://ptrs-ukb/tmp-output/hail_prs_test_quick.log \
    >> $JOBNAME.log 2>&1 

# --gwas-ht gs://ptrs-ukb/results/gwas_on_subsets/all_chrs_up.ht \
# hailctl dataproc stop $JOBNAME

