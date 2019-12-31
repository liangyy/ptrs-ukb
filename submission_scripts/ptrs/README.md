Run on CRI.

Example:

```
indivlist=British-test
for i in `seq 1 17`
do 
  bash submit.sh $i $indivlist-$i ctimp
done
```

To specify prediction model:

```
indivlist=British-test
for i in `seq 1 17`
do 
  bash submit.sh $i $indivlist-$i mesa AFA
done
```
