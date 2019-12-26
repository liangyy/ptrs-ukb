Run on CRI.

Example:

```
indivlist=British-test
for i in `seq 1 17`
do 
  echo $i
  bash submit.sh $i $indivlist-$i
done
```
