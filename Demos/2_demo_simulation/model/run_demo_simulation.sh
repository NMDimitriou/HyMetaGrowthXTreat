#! /bin/bash -f

#nvidia-smi --query-gpu=index,timestamp,name,driver_version,pstate,utilization.gpu,utilization.memory,memory.total,memory.used --format=csv,nounits --filename=gpu-log.csv
#sed -i 's/,\s/\t/g' gpu-log.csv
#less -S gpu-log.csv | awk -F'\t' '$9 < 16000 {print $9}' > gpu_id.txt


echo "------------------------------------------------------"
echo "          RUNNING DATASET DEMO                        "
echo "------------------------------------------------------"
cd demo_simulation/model/parent/
./main_parent start ../../../IC/corr_dens corr_dens
cd ../
./main ../params.txt corr_dens 1
cd parent/
./main_parent stop ../../../IC/corr_dens corr_dens
cd ../
make clear


