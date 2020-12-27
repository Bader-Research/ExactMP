#include "header.h"

int main(int argc, char **argv) {

   char *opFile;

   time_t st0, st1, et0, et1, st2,
          ut0, ut1, ft0, ft1, st3;

   checkUsage(argc, argv);

   opFile = (char *)calloc(20, sizeof(char));
   strcat(opFile, argv[1]);
   strcat(opFile, "Vad");
   fp = fopen(opFile,"w");
   assert(fp);
   
   st0 = time((time_t *)NULL);
   st1 = st0;
   st2 = st0;
   st3 = st0;

   printf("\n\t\tSTARTING JOB %s @ TIME: %s",argv[1], ctime(&st2));
   printf("\t     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
   fprintf(fp,"\n\t\tSTARTING JOB %s @ TIME: %s",argv[1], ctime(&st3));
   fprintf(fp,"\t     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");

   printf("\nSetting framework...\n");
   setFrameWork(argv[1]);
   
   ut0 = time((time_t *)NULL);
   ut1 = ut0;

   printf("\nPNI Sites: %d, Const Sites: %d, PI Sites: %d\n",numPNISites, numConstSites, matrix.num_pars_inf_sites);
   fprintf(fp,"\nPNI Sites: %d, Const Sites: %d, PI Sites: %d\n",numPNISites, numConstSites, matrix.num_pars_inf_sites);
   printf("\nApplying B&B with a score of %d...\n",bestCost+pNiCost);
   fprintf(fp,"\nApplying B&B with a score of %d...\n",bestCost+pNiCost);

   ft0 = time((time_t *)NULL);
   ft1 = ft0;
   
   if (matrix.num_taxa < 8) {
      branchAndBound();
   } else {
      bnb();
   }

   et0 = time((time_t *)NULL);
   et1 = et0;

   printf("\n         ^^^^^^^^^^^^^^^^^^^^^\n");
   printf("\t   RESULT STATISTICS");
   printf("\n         ^^^^^^^^^^^^^^^^^^^^^\n");
   printf("1.Number of trees retained with best scores: %d\n",numBestSolutions);
   printf("2.B&B started with an initial score of %d\n",initBestCost+pNiCost);
   printf("3.B&B converged to the best score of %d\n",bestCost+pNiCost);
   printf("4.Start/End time for setting up B&B framework\n");
   printf("    Start  time: %s",ctime(&st0));
   printf("    Finish time: %s",ctime(&ut0));
   printf("5.Start/End time of actual B&B algorithm\n");
   printf("    Start  time: %s",ctime(&ft0));
   printf("    Finish time: %s",ctime(&et0));
   printf("\n\t\t-o-o-o-\n\n");

   fprintf(fp,"\n         ^^^^^^^^^^^^^^^^^^^^^\n");
   fprintf(fp,"\t   RESULT STATISTICS");
   fprintf(fp,"\n         ^^^^^^^^^^^^^^^^^^^^^\n");
   fprintf(fp,"1.Number of trees retained with best scores: %d\n",numBestSolutions);
   fprintf(fp,"2.B&B started with an initial score of %d\n",initBestCost+pNiCost);
   fprintf(fp,"3.B&B converged to the best score of %d\n",bestCost+pNiCost);
   fprintf(fp,"4.Start/End time for setting up B&B framework\n");
   fprintf(fp,"    Start  time: %s",ctime(&st1));
   fprintf(fp,"    Finish time: %s",ctime(&ut1));
   fprintf(fp,"5.Start/End time of actual B&B algorithm\n");
   fprintf(fp,"    Start  time: %s",ctime(&ft1));
   fprintf(fp,"    Finish time: %s",ctime(&et1));
   fprintf(fp,"\n\t\t-o-o-o-\n\n");

   fclose(fp);
   
   printf("\nClearing memory...\n");
   if (matrix.num_taxa < 8) {
      clearMemory();
   } else {
      //clearWorkSpace();
   }

   return 0;
}

void checkUsage(int argc, char **argv) {
   if (argc!=3) {
     printf("USAGE: vc {seq_matrix_file} {arrangements_per_level}\n");
     exit(0);
   }
   randOptLev = 1;
   NUM_TBR_TREES = atoi(argv[2]);
}
