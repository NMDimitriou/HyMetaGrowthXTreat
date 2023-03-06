#!/bin/bash
#SBATCH --account=def-gmitsis
#SBATCH --job-name=CAL_KSC_ALEE_PAC_AC_EXP_IC_SET_2_RUN_1_NO_R
#SBATCH --output=%x-%j.out
#SBATCH --time=0-20:00
#SBATCH --mail-user=nikolaos.dimitriou@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --gpus-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH -p compute_full_node

cd $SLURM_SUBMIT_DIR
module load cuda/11.0.3 gcc/9.4.0 gsl/2.7 spectrum-mpi/10.4.0
scontrol show hostnames > nodelist-$SLURM_JOB_ID
export PATH=$HOME/usr/torc/bin:$PATH

echo "------------------------------------------------------"
echo "          RUNNING DATASET AC SERIES 2            	    "
echo "------------------------------------------------------"
cd model/parent/
mpirun --map-by ppr:1:node ./main_parent start ../../IC/series2_treatment/Pac_0p5_AC Pac_0p5_AC
#./main_parent start.txt avg_s1
cd ../../

sleep 5

mpirun -hostfile nodelist-$SLURM_JOB_ID ./sample
#mpirun -np 64 ./sample
#dakota sa_lhs.in

sleep 5

cd parent/
mpirun --map-by ppr:1:node ./main_parent stop ../../IC/series2_treatment/Pac_0p5_AC Pac_0p5_AC
#./main_parent stop.txt avg_s1
cd ../../

mkdir $SLURM_JOB_NAME
#mv dakota* *.out runs $SLURM_JOB_NAME
mv *.txt *.out $SLURM_JOB_NAME
