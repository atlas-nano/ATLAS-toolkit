"""
#--------------------------
# Name: setup_cpu_bench.py
# Purpose: Automate the creation of slurm files and directory structure for LAMMPS simulations benchmarking.
# 
# The script creates a "benchmarking" directory with subdirectories named after the number of tasks per slurm job.
# It copies the lammps input and data files into each subdirectory.
# For each directory, it also creates a slurm file with the necessary parameters filled.
# 
# Input Parameters:
# input_file, REQUIRED(str): The input file for the simulation with the format 'in.prefix'.
# time (str): The time limit for each slurm job, in the format 'hh:mm:ss'.
# account (str): The account name to be billed for resource usage.
# rtemp (float): The temperature for the simulation in Kelvins.
# press (float): The pressure for the simulation in atmospheres.
#
# This script assumes that the lammps input and data files are present in the working directory where it's run from, and the respective names are 'in.prefix' and 'data.prefix'.
#--------------------------
"""
import os
import shutil
import argparse

# Define directory from which script is run

# List of tasks numbers
def setup_benchmark(input_file, account, ntasks_values, time, rtemp, press):
    # Create the benchmarking directory

    main_dir = os.getcwd()

    os.makedirs("benchmarking", exist_ok=True)

    prefix = input_file.split(".")[1]

    for ntasks in ntasks_values:
        # Create subdirectory named after the number of tasks
        subdirectory_path = os.path.join("benchmarking", str(ntasks) + "_tasks")
        os.makedirs(subdirectory_path, exist_ok=True)

        # Copy lammps data and input files
        shutil.copy2(os.path.join(main_dir, 'in.' + prefix), os.path.join(subdirectory_path, 'in.' + prefix))
        shutil.copy2(os.path.join(main_dir, 'data.' + prefix), os.path.join(subdirectory_path, 'data.' + prefix))

        # Create slurm file
        slurm_file_path = os.path.join(subdirectory_path, "slurm_" + str(ntasks) + ".sh")

        slurm_text = f"""#!/bin/bash
#SBATCH --job-name={prefix}_job
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={ntasks}
#SBATCH --time={time}
#SBATCH --output={prefix}_%j.out
#SBATCH --account={account}
# Load necessary modules
module purge
module load slurm
module load cpu/0.15.4
module load gcc/9.2.0
module load openmpi
module load cmake
module load openblas
module load amdfftw
module load gsl
module load netlib-scalapack
module load netlib-lapack
# File names
lmp_equil_file=in.{prefix}
lmp_data_file=data.{prefix}
lmp_log_file={prefix}.log
# LAMMPS executable and parameters
PARALLEL="mpirun -n {ntasks}"
LMP="/home/rramji/codes/LAMMPS/speed_lammps/build_custom/lmp_rob -var rtemp {rtemp} -var press {press}"
# Echo job details
echo "LAMMPS dynamics of {prefix} at {rtemp}K"
# Run LAMMPS
$PARALLEL $LMP -in ${{lmp_equil_file}} -log ${{lmp_log_file}} -var rtemp {rtemp} -var press {press}
"""
        with open(slurm_file_path, 'w') as f:
            f.write(slurm_text)

    print("Benchmarking directories and files were created successfully.")

def create_batch_submission_script(ntasks_values):
    bash_submit_script_path = os.path.join("benchmarking", "submit_all.sh")

    with open(bash_submit_script_path, 'w') as f:
        for ntasks in ntasks_values:
            subdirectory_path = os.path.join(str(ntasks) + "_tasks")
            slurm_file_name = "slurm_" + str(ntasks) + ".sh"
            f.write(f"cd {os.path.join(subdirectory_path)}\n")
            f.write(f"sbatch {os.path.join(slurm_file_name)}\n")
            f.write(f"cd ..\n")

    print(f"Bash script to submit all slurm jobs created at {bash_submit_script_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Set up subdirectories for LAMMPS benchmarking simulations")
    parser.add_argument("input_file", type=str, help="LAMMPS input file")
    parser.add_argument("--account", type=str, default="csd799", help="project/account info (eg. csdXXX or ddpXXX)")
    parser.add_argument("--ntasks_values", type=str, default=[4, 8, 16, 32, 64], help="Number of tasks (MPI processes)")
    parser.add_argument("--time", type=str, default="08:00:00", help="Run time (hh:mm:ss)")
    parser.add_argument("--rtemp", type=int, default=300, help="Temperature")
    parser.add_argument("--press", type=int, default=1, help="Pressure")


    args = parser.parse_args()

    setup_benchmark(args.input_file, args.account, args.ntasks_values, args.time, args.rtemp, args.press)
    create_batch_submission_script(args.ntasks_values)
