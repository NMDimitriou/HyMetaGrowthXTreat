#!/bin/bash
#SBATCH --account=def-gmitsis
#SBATCH --job-name=run2_TREAT_HYBRID_NO_R_adhes-5_phenchange0
#SBATCH --output=%x-%j.out
#SBATCH --time=0-04:30
#SBATCH --mail-user=nikolaos.dimitriou@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1

cd $SLURM_SUBMIT_DIR
#module load nixpkgs/16.09 gcc/7.3.0  gsl/2.5 cuda/10.1 openmpi/3.1.2
module load gcc cuda

make

s2=(Pac_0p5_AC Pac_0p5_AN Pac_0p5_AS Pac_0p5_AW Pac_0p5_AE Pac_0p05_BC Pac_0p05_BN Pac_0p05_BS Pac_0p05_BW Pac_0p05_BE Pac_0p005_CC Pac_0p005_CN Pac_0p005_CS Pac_0p005_CW Pac_0p005_CE Pac_0p0005_DC Pac_0p0005_DN Pac_0p0005_DS Pac_0p0005_DW Pac_0p0005_DE)

#s2=(Pac_0p005_CN Pac_0p005_CS Pac_0p005_CW Pac_0p005_CE Pac_0p0005_DC Pac_0p0005_DN Pac_0p0005_DS Pac_0p0005_DW Pac_0p0005_DE)

s1=(AC AN AS AW AE BC BN BS BW BE CC CN CS CW CE DC DN DS DW DE)
#s1=(CN CS CW CE DC DN DS DW DE)
j=0

for i in ${s2[@]}
do
    echo "------------------------------------------------------"
    echo "          RUNNING DATASET s2 $i	    	            "
    echo "------------------------------------------------------"
    ./main best_params_treatment/final_${s1[$j]}s2_best.txt ${i} 1
	((j++))
done

mkdir -p $SLURM_JOB_NAME
mv *.txt *.csv *.out $SLURM_JOB_NAME
