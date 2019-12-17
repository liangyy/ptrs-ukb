## About

Different from `scripts/`, here we keep submission scripts using pipelines built at `pipeline`. 
The main difference is that the jobs here contain multiple steps.

And here contains the submission scripts for GCP.

## How do you mean by submission?

This repo is mainly used for running jobs on `nucleus` where no cluster scheduler is in use. 
Instead the jobs is submitted by using `screen`.

To make it resemble using a cluster scheduler and simplify the work, please stick with the following rules for 'submitting' a job.

First of all, the job should be written in a script just like the job script you would write for PBS or Slurm.
Here the scripts named `*.screen` should follow the rules strictly.
This means that you need something as follow

```
# set up computing environment
# conda activate ...
# export ...

# specify the working directory
# cd ...

# content of your job 
# cmd1 
# cmd2 
# cmd3 > log

# at the end of the job, you need to add this extra line to kill the screen
screen -X kill
```

**To submit the job**: 

Open a new terminal and

```
$ screen -dmS job_name bash script_name.screen
```

**To monitor screen**:

```
$ screen -ls
```
