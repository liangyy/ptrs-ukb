To submit, an example run is

```
outdir=/vol/bmd/yanyul/UKB/predicted_expression
for p in `cat pop-list-minimal.txt`
do
  screen -dmS $p bash format_pred_expr.screen $p $outdir
done
```
