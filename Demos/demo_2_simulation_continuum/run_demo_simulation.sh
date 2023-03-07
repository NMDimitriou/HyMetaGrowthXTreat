#! /bin/bash -f

echo "------------------------------------------------------"
echo "          RUNNING DATASET DEMO                        "
echo "------------------------------------------------------"
cd parent/
./main_parent start ../../IC/ new_dens_D5
cd ../
./main params_treatment.txt new_dens_D5 1
cd parent/
./main_parent stop  ../../IC/ new_dens_D5
cd ../
make clear
./kill_ipcs.sh
