LD clumping using plink1.9

```
plink \
  --bgen {input[0]} \
  --clump {input[1]} \
  --clump-p1 1 \
  --clump-r2 0.1 \
  --clump-kb 250 \
  --clump-snp-field {config[gwas][snp-field]} \
  --clump-field {config[gwas][pval-field]} \
  {memory_cmd} \
  {thread_cmd} \
  --out {wildcards.outdir}/ld-clumping.{wildcards.name}.chr{wildcards.i}
```

Here we have 289 GWASs to run LD-clumping.
Each GWAS clumping job is split into 22 sub-jobs.
And they get to merged as one at the end.

Unfortunately, the bgen files from `subset_bgen.snmk` are not directly usable for LD-clumping due to two issues:

1. bgen v1.2 is not readable for plink1.9
2. there are duplicated variant ID (rsID) 

So, we need to run LD-clumping in two steps:

1. run `cleanup.snmk` to convert BGEN v1.2 to VCF without duplicated variant
2. run `ld_clump.snmk`
