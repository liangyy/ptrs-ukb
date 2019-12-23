if [[ ! -d logs ]]
then
  mkdir logs
fi

for i in `cat ../../misc/trait_list.txt`
do
  qsub -v TRAIT=$i -N spredixcan-$i run.qsub
done