# ARGS1: subset number
# ARGS2: individual list 

NUM=$1
INDIVLIST=$2

if [[ ! -d logs ]]
then
  mkdir logs
fi


spredixcanList=`for i in `cat external_data/martin_et_al_2019ng_table_s6.csv |cut -f 1 -d,|awk '{print tolower($0)}'|tail -n +2`; do echo 'subset$NUM_x_'$i; done|tr '\n' ','`

echo qsub -v INDIVLIST=$INDIVLIST,SPREDIXCANLIST=$spredixcanList,NUM=$NUM -N $NUM-$INDIVLIST run.qsub
