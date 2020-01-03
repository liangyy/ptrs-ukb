# ARGS1: CONFIGFILE middle name
# ARGS2: MODELTAG 

if [[ ! -d logs ]]
then
  mkdir logs
fi

for i in `cat ../../misc/trait_list.txt`
do
  qsub -v TRAIT=$i -N annot-$i run_annot.qsub
done
