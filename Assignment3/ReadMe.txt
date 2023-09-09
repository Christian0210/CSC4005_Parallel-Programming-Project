The command line of execution for each file:

hw3_seq.cpp                       g++ hw3_seq.cpp -o seq.out -lX11 && ./seq.out
 
hw3_mpi.cpp                      mpic++ hw3_mpi.cpp -o mpi.out -lX11 && mpirun -n 4 ./mpi.out

hw3_pthread.cpp                g++ hw3_pthread.cpp -o pthread.out -lX11 -lpthread && ./pthread.out 4

hw3_openmp.cpp                g++ hw3_openmp.cpp -o openmp.out -lX11 -fopenmp && ./openmp.out 4

hw3_mpi_openmp.cpp        mpic++ hw3_mpi_openmp.cpp -o mpi_openmp.out -lX11 -fopenmp && mpirun -n 4 ./mpi_openmp.out 4


