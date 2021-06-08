traits=`cat trait_list.txt`

for trait in $traits
do
  echo qsub -v TRAIT=$trait -N run_prs_$trait run_prs_by_trait.qsub
done
