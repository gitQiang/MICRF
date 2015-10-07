#! /bin/bash
#$ -l mem=9G,time=:30: -cwd -S /bin/bash -N test

#matlab -nodisplay -r "qsub($SGE_TASK_ID);exit;"

onefile="iRef_max_1.mat"
dir="iRefbest/"
matlab -nodisplay -r "qsub($SGE_TASK_ID,'$onefile','$dir');exit;"

#matlab -nodisplay -r "qsubtest($SGE_TASK_ID,'$onefile','$dir');exit;"
