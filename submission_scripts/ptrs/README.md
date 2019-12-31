Run on CRI.

# Compute PTRS

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
model=AFA
for i in `seq 1 17`
do 
  bash submit.sh $i $indivlist-$i mesa $model
done
```

For other population, it is recommended to run in two steps

* Step 1

```
model=AFA
for i in `echo "African Chinese Indian"`
do
  bash submit.sh 1 $i mesa $model
done
```

* Step 2

```
model=AFA
indivlist=African
for i in `seq 2 17`
do
  bash submit.sh $i $indivlist mesa $model
done
```

# Compute R2

Example:

```
bash submit_r2.sh r2
```

Run MESA models

```
bash submit_r2_mesa.sh config.r2_mesa_template.yaml CAU
```
