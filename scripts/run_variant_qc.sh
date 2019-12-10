export PYSPARK_SUBMIT_ARGS="--driver-memory 50G --executor-memory 10g pyspark-shell
conda activate /vol/bmd/yanyul/softwares/conda_envs/hail

WORKDIR=/homes/yanyul/labshare/GitHub/ptrs-ukb
OUTPREF=/vol/bmd/yanyul/UKB/variant_qc/imp_all
cd $WORKDIR


python code/variant_qc.py \
  --bgen-path /vol/bmd/data/ukbiobank/genotypes/v3/ukb_imp_chr{chr_num}_v3.bgen \
  --bgen-sample /vol/bmd/data/ukbiobank/genotypes/v3/ukb19526_imp_chr1_v3_s487395.sample \
  --bgen-index /vol/bmd/yanyul/UKB/bgen_idx/ukb_imp_chr{chr_num}_v3.bgen.idx2 \
  --hail-block-size 512 \
  --output-filename $OUTPREF \
  > $OUTPREF.log 2>&1
  
screen -X kill