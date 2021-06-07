# args1: trait file
# args2: bgen prefix
# args3: range file
# args4: out prefix
# args5: extra opts

plink2=/gpfs/data/im-lab/nas40t2/yanyul/softwares/plink2

$plink2 --bgen $2.bgen --sample $2.sample \
  --score $1 1 3 4 header-read cols=scoresums \
  --q-score-range $3 $1 header \
  --out $4 $5
