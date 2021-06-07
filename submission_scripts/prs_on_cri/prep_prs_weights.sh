outdir=/gpfs/data/im-lab/nas40t2/yanyul/PTRS/prs_weights

mkdir $outdir

Rscript ../../scripts/gen_prs_weight_table.R \
  --data_yaml prs_weight.yaml \
  --pval_cutoffs 5e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,0.05,0.1,0.5,1 \
  --snpid rsid \
  --pval pval \
  --effect_allele alt \
  --effect_size beta \
  --output_prefix $outdir/ld_clump. > prep_prs_weights.log 2>&1
  