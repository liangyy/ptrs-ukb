Example run:

```
model=ctimp_Whole_Blood
outdir=/vol/bmd/yanyul/UKB/predicted_expression
screen -dmS $model bash run_pred_expr.screen $model $outdir
```

Run MESA models:

```
model=AFA  # AFHI CAU ALL HIS
outdir=/vol/bmd/yanyul/UKB/predicted_expression
screen -dmS $model bash run_pred_expr_mesa.screen $model $outdir
```

