#PBS -l walltime=00:45:00
#PBS -l select=4:ncpus=8:mpiprocs=8:mem=64gb
#PBS -N NEGFjob
#PBS -o job.out
#PBS -e job.err
#PBS -A IscrC_GaNLEDs

cd $PBS_O_WORKDIR
module load autoload intelmpi/5.0.2--binary
module load autoload boost/1.58.0--intelmpi--5.0.2--binary
module load autoload mkl/11.2--binary
module load autoload fftw/3.3.4--intelmpi--5.0.2--binary

# this one for intelmpi
mpirun -np 32 ./angel.bin pn ../results/try.log
