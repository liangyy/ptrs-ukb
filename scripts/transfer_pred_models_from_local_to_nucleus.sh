WORKDIR=/Users/yanyul/Documents/repo/github/ptrs-ukb
SOURCEDIR=/Users/yanyul/Desktop/tmp/

rsync -avc --log-file=$SOURCEDIR/logs/to-nucleus.log $SOURCEDIR/gtex_v8_pred_models/ nucleus:/vol/bmd/yanyul/data/gtex_v8_pred_models

screen -X kill
