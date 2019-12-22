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
