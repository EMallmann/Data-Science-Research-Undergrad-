#!/bin/bash

#SBATCH --partition=batch             #Partition to submit to
#SBATCH --time=6-00:00:00             #Time limit for this job
#SBATCH --nodes=1                     #Nodes to be used for this job during runtime. This shouldn't ever be more than one unless one can prove it is necessary.
#SBATCH --ntasks-per-node=1           #Number of CPUs. Cannot be greater than number of CPUs on the node. This shouldn't ever be more than four unless one can prove it is necessary.
#SBATCH --nodelist=compute67                 #for example, since 68 is idle at the moment
#SBATCH --mem=70GB                     #Total memory for this job
#SBATCH --job-name="WorkChr1_23 RandomForest Computing"     #Name of this job in work queue
#SBATCH --output=chr1_23_work_randomForest-ssample.out          #Output file name
#SBATCH --error=chr1_23_work_randomForest-ssample.err          #Error file name
#SBATCH --mail-user=foongmin@hotmail.com  #Email to send notifications to
#SBATCH --mail-type=ALL               #Email notification type (BEGIN, END, FAIL, ALL). To have multiple use a comma separated list. i.e END,FAIL.


ulimit -s 100000
module load R
Rscript --max-ppsize=500000 chr1_23_randomForest_work.R

#Rscript --max-mem-size=500M --max-ppsize=5000000 --max-vsize-500M chr1_23_randomForest.R
