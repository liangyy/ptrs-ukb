NAME=test1
hailctl dataproc start $NAME-chr22 -p 40 --preemptible-worker-boot-disk-size 16 > $NAME\_chr22.log 2>&1 

hailctl dataproc submit $NAME-chr22 --properties spark.driver.memory=64g,spark.executor.memory=10g,spark.master=cluster,spark-deploy= --pyfiles ../../code/my_hail_helper.py,../../code/gwas_helper.py ../../code/gwas_on_subset_ht.py --bgen-path gs://ukb_data/genotypes/v3/ukb_imp_chr{chr_num}_v3.bgen --bgen-sample gs://ukb_data/genotypes/ukb19526_imp_chr1_v3_s487395.sample --bgen-index gs://ukb_data/bgen_idx/ukb_imp_chr{chr_num}_v3.bgen.idx2 --hail-block-size 512 --chrs 22 --output-filename gs://ptrs-ukb/tmp-output/test_run/$NAME\_chr22  --variant-ht gs://ukb_data/variant_qc/imp_all.ht --variant-qc-yaml gs://ptrs-ukb/results/misc/gwas_full_variant_qc.yaml --pheno-covar-path gs://ptrs-ukb/results/pheno_and_covar/all_subsets.{subset_name}.ht --subset-yaml gs://ptrs-ukb/results/misc/gwas_full_subset_list.yaml --google-cloud-project ukb-im >> $NAME\_chr22.log 2>&1 

hailctl dataproc stop $NAME-chr22

screen -X kill