#! make some required directories
# mkdir -p ENV
mkdir -p RESULTS
mkdir -p PLOTS

#! remove the old and useless files
find . -name "*.exe" -type f -delete
find . -name "*.x" -type f -delete
find . -name "*.txt" -type f -delete

#! run serial program 
g++ LLF.cpp -o prg_serial.x
./prg_serial.x
python3 plot_serial.py

#! remove the old and useless files
find . -name "*.exe" -type f -delete
find . -name "*.x" -type f -delete
find . -name "*.txt" -type f -delete

#! run parallel program with 2 cores
mpicxx LLF_mpi.cpp -o prg_mpi_2.x
mpirun -N 2 ./prg_mpi_2.x
python3 plot_parallel.py

#! remove the old and useless files
find . -name "*.exe" -type f -delete
find . -name "*.x" -type f -delete
find . -name "*.txt" -type f -delete

#! run parallel program with 2 cores
mpicxx LLF_mpi.cpp -o prg_mpi_4.x
mpirun -N 4 ./prg_mpi_4.x
python3 plot_parallel.py

