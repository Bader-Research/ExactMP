#include "header.h"

int main(int argc, char **argv) {

   char *opFile, *inFile;
   double startTime, endTime;
   int ctr, finNum;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

   startTime = MPI_Wtime();

   if (my_rank==0) {
     printf("argc: %d\n",argc);
     fflush(stdout);
     for (ctr=0; ctr < argc; ctr++) {
      printf("argv[%d]: %s\n",ctr,argv[ctr]);
      fflush(stdout);
     }
   }
   checkUsage(argc, argv);


   //printf("\nSetting framework...\n");
   inFile = (char *)calloc(200,sizeof(char));
   strcat(inFile,"/tmp/");
   strcat(inFile,argv[1]);
   setFrameWork(inFile);

   if (my_rank==0) {
     printf("\nRank: %d, PNI Sites: %d, Const Sites: %d, PI Sites: %d\n",my_rank, numPNISites, numConstSites, matrix.num_pars_inf_sites);
     fflush(stdout);
   }
   if (my_rank==0) {
     printf("\nRank: %d, Applying B&B with a score of %d...\n",my_rank, bestCost+pNiCost);
     fflush(stdout);
   }

   //bestCost = 2000053612;
   if (matrix.num_taxa < 8) {
      branchAndBound();
   } else {
      //bnb();
      paBnB();
   }

   /*
   if (my_rank==0) {
     printf("Rank: %d, B&B converged to the best score of %d\n",my_rank, bestCost+pNiCost);
     fflush(stdout);
   } */

/*
   printf("\nClearing memory...\n");
   if (matrix.num_taxa < 8) {
      clearMemory();
   } else {
      clearWorkSpace();
   }
*/

   printf("Rank: %d, B&B converged to the best score of %d\n",my_rank, bestCost+pNiCost);
   fflush(stdout);

   endTime = MPI_Wtime();

   printf("Time: %.2f\n",endTime-startTime);

   MPI_Finalize();
   //printf("Rank: %d, Finalize: %d\n",my_rank,finNum);

   return 0;
}

void checkUsage(int argc, char **argv) {
   randOptLev = 3;
   if (argc==3) {
     NUM_TBR_TREES = atoi(argv[2]);
   } else {
     printf("USAGE: ExactMP {seq_matrix_file} {arrangements_per_level}\n");
     exit(0);
   }
}
