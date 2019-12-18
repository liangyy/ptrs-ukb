## About

Here we run REML-based regulability using two implementations: `gcta` and `hail`.

1. First, we run `run_format.screen` which performs formatting; 
2. Then, we run `run_gcta.screen` which runs `gcta` based approach;
3. In parallel to 2, we run `run_hail.screen` which runs `hail` based approach.

## Example run

For a given individual list, set it up.

```
POPNAME='British-test-1'  
OUTDIR=/vol/bmd/yanyul/UKB/gcta_regulability
```

Do formatting first.

```
screen -dmS format-$POPNAME bash run_format.screen $POPNAME $OUTDIR
```

Then, do gcta

```
screen -dmS format-$POPNAME bash run_gcta.screen $POPNAME $OUTDIR
```