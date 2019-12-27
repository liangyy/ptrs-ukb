# ARGS1: subset number
# ARGS2: individual list 

NUM=$1
INDIVLIST=$2

if [[ ! -d logs ]]
then
  mkdir logs
fi

traitlist=`cat ../../external_data/martin_et_al_2019ng_table_s6.csv |cut -f 1 -d,|awk '{print tolower($0)}'|tail -n +2`
spredixcanList=`for i in $traitlist; do echo 'subset'$NUM'_x_'$i; done|tr '\n' ':'|sed 's#:$##g'`

qsub -v INDIVLIST=$INDIVLIST,SPREDIXCANLIST=$spredixcanList,NUM=$NUM -N $NUM-$INDIVLIST run.qsub
