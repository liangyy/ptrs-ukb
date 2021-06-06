# args1: bgen
# args2: sample
# args3: outprefix
# args4: RAM in Mb
# args5: threads

plink2=/gpfs/data/im-lab/nas40t2/yanyul/softwares/plink2

$plink2 --bgen $1 --sample $2 --make-grm-bin --out $3 --memory $4 --threads $5
