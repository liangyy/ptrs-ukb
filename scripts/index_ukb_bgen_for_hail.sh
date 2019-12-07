# # environment setups:
# export PYSPARK_SUBMIT_ARGS="--driver-memory 40G --executor-memory 8g pyspark-shell
# conda activate /vol/bmd/yanyul/softwares/conda_envs/hail

WORKDIR=/vol/bmd/yanyul/UKB/bgen_idx
BGENFOLDER=/vol/bmd/data/ukbiobank/genotypes/v3
IDXFOLDER=/vol/bmd/yanyul/UKB/bgen_idx
CODEFOLDER=/vol/bmd/yanyul/GitHub/ptrs-ukb/code
SAMPLEFILE=ukb19526_imp_chr1_v3_s487395.sample
LOGFILE=index_ukb_bgen_for_hail.log

for i in `seq 1 22`
do 
  python $CODEFOLDER/index_bgen_for_hail.py \
    --bgen ukb_imp_chr$i\_v3.bgen \
    --index_file ukb_imp_chr$i\_v3.bgen.idx2 \
    >> $WORKDIR/$LOGFILE 2>&1
done
