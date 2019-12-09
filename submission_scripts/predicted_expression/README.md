To submit, an example run is

```
outdir=/vol/bmd/yanyul/UKB/predicted_expression
for p in `ls pop-list-minimal.txt`
do
  screen -S $p bash format_pred_expr.screen $p $outdir
done
```