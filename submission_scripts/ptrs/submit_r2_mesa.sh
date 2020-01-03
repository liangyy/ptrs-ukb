# ARGS1: configfile template file path
# ARGS2: MESA model name

if [[ ! -d logs ]]
then
  mkdir logs
fi

TEMPLATEPATH=$1
MESANAME=$2

OUT=config.mesa-$MESANAME.yaml

cat $TEMPLATEPATH | sed "s#MESA-MODEL-NAME#$MESANAME#g" > $OUT

for i in `seq 1 17`
do
  qsub -v NUM=$i,CONFIG=mesa-$MESANAME -N R2-$i-$MESANAME run_r2.qsub
done
