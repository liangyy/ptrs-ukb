WORKDIR=/Users/yanyul/Documents/repo/github/ptrs-ukb
TARGETDIR=/Users/yanyul/Desktop/tmp
# cd $WORKDIR

mkdir -p $TARGETDIR/gtex_v8_pred_models
mkdir -p $TARGETDIR/logs

mylist=`cat $WORKDIR/misc/list_of_gtex_v8_pred_models.txt`


for i in $mylist
do 
  names=`echo $i | awk '{n=split($1,a,"/"); print a[n-1]"/"a[n]}'`
  tag=`echo $names | tr '/' '_'`
  mkdir -p $TARGETDIR/gtex_v8_pred_models/$names
  rsync -avc --log-file=$TARGETDIR/logs/$tag.log gardner:/gpfs/data/im-lab/nas40t2/abarbeira/projects/gtex_v8/models_v1/$names/ $TARGETDIR/gtex_v8_pred_models/$names
done

screen -X kill
