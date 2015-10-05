#! /bin/bash
#$ -l mem=5G,time=:20: -cwd -S /bin/bash -N test

matlab -nodisplay -r "qsub($SGE_TASK_ID);exit;"
