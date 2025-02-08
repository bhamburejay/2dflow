/*The Parallel Hello World Program*/
#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv)
{
   int node, size;
   
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &node);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
     
   printf("Hello World from Node %d\n",node);
            
   MPI_Finalize();
   return 0;
}