source /vol/bmd/yanyul/miniconda3/etc/profile.d/conda.sh
conda activate /vol/bmd/yanyul/softwares/conda_envs/hail

python ../code/ld_clump_prs.py \
  --bgen-path /vol/bmd/data/ukbiobank/genotypes/v3/ukb_imp_chr{chr_num}_v3.bgen \
  --bgen-sample /vol/bmd/data/ukbiobank/genotypes/v3/ukb19526_imp_chr1_v3_s487395.sample \
  --bgen-index /vol/bmd/yanyul/UKB/bgen_idx/ukb_imp_chr{chr_num}_v3.bgen.idx2 \
  --chrs 22 \
  --subset-and-gwas ../misc/prs_test.yaml \
  --variant-pool /vol/bmd/yanyul/UKB/gwas_on_subset/gwas_runs_in_tsv/gwas_in_tsv_subset3_x_height.tsv \
  --pval-thresholds 5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,0.05,0.1,0.5,1 \
  --output-prefix /vol/bmd/yanyul/UKB/tmp/test_prs \
  --hail-log /vol/bmd/yanyul/UKB/tmp/test_prs.log
  