# ARGS1: CONFIGFILE middle name
# ARGS2: MODELTAG 

if [[ ! -d logs ]]
then
  mkdir logs
fi

for i in `cat ../../misc/trait_list.txt`
do
  qsub -v TRAIT=$i,CONFIG=$1,MODELTAG=$2 -N spredixcan-$i run.qsub
done
