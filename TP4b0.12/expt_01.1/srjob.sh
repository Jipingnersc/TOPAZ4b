#!/bin/bash -l

#SBATCH --account=nn2993k

#SBATCH -J  TP4_B2

#SBATCH -N 7   # number of nodes
#SBATCH -n 224   # number of cores 
## Each compute node has 24 cores (See more details in section Hardware on Sisu User Guide).
## Calculate the total number of cores and store it in variable ncores
##(( ncores = SLURM_NNODES * 24 ))
##(( ncores = SLURM_NNODES * 32 ))  #on Fram


##SBATCH -p small_long
## Choose a suitable queue <test,small,large>
## How to check queue limits: scontrol show part <queue name>
## for example: scontrol show part small

## System message output file
#SBATCH -o log/HYCICE.%J.out

## System error message file
#SBATCH -e log/HYCICE.%J.err

## How long job takes, wallclock time hh:mm:ss
#SBATCH -t 00:18:00
##SBATCH -t 70:00:00

##SBATCH  --mail-type=END
##SBATCH --mail-user=jiping.xie@nersc.no

## option: -n (total number of mpi processes)
## option: -N (number of mpi processes per compute node)
## option: -S (number of mpi processes per NUMA node)
## option: -ss (allocate memory only from a local NUMA node)


export SLURM_SUBMIT_DIR=$(pwd)
# Enter directory from where the job was submitted
cd $SLURM_O_WORKDIR       ||  { echo "Could not go to dir $SLURM_O_WORKDIR  "; exit 1; }

# Enter directory from where the job was submitted
cd $SLURM_SUBMIT_DIR       ||  { echo "Could not go to dir $SLURM_O_WORKDIR  "; exit 1; }


# Initialize environment (sets Scratch dir ($S), Data dir $D ++ )
source ../REGION.src  || { echo "Could not source ../REGION.src "; exit 1; }
source ./EXPT.src  || { echo "Could not source EXPT.src"; exit 1; }
echo "NMPI =$NMPI (Number of MPI tasks needed for running job) "

### Enter Scratch/run dir and Run model
cd $S || { echo "Could not go to dir $S  "; exit 1; }

    #aprun -n $NMPI -m 1000M ./hycom ${e} # Run hycom
srun --mpi=pmi2 -n $NMPI --cpu_bind=cores ./hycom 

    # (A) actually, we want the batch to continue even if one of the runs fails
    #
    # if [ $? -ne 0 ]
    # then
    #     exit $?
    # fi

