# ARGS: if you need real run, add something here
#       otherwise, it will be dry run

if [[ ! -z $1 ]]
then 
  rsync -avc /Users/yanyul/Documents/repo/github/ptrs-ukb/output/ nucleus:/vol/bmd/yanyul/GitHub/ptrs-ukb/output
else
  rsync -navc /Users/yanyul/Documents/repo/github/ptrs-ukb/output/ nucleus:/vol/bmd/yanyul/GitHub/ptrs-ukb/output
fi
