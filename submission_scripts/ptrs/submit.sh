# ARGS1: subset number
# ARGS2: individual list 
# ARGS3: configfile middle name
# ARGS4: (optional) predicted expression model 

NUM=$1
INDIVLIST=$2
CONFIG=$3
if [[ ! -z $4 ]]
then
  PREDMODEL=$4
fi

if [[ ! -d logs ]]
then
  mkdir logs
fi

traitlist=`cat ../../external_data/martin_et_al_2019ng_table_s6.csv |cut -f 1 -d,|awk '{print tolower($0)}'|tail -n +2`
spredixcanList=`for i in $traitlist; do echo 'subset'$NUM'_x_'$i; done|tr '\n' ':'|sed 's#:$##g'`

if [[ -z $PREDMODEL ]]
then
  qsub -v INDIVLIST=$INDIVLIST,SPREDIXCANLIST=$spredixcanList,NUM=$NUM,CONFIG=$CONFIG -N $NUM-$INDIVLIST-$CONFIG run.qsub
else
  qsub -v INDIVLIST=$INDIVLIST,SPREDIXCANLIST=$spredixcanList,NUM=$NUM,CONFIG=$CONFIG,PREDMODEL=$PREDMODEL -N $NUM-$INDIVLIST-$CONFIG-$PREDMODEL run.qsub
fi
