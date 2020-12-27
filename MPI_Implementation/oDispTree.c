#include "header.h"

void displayTree(int **cst) {
  int ctr=0, tmp=0;

  printf("\n\t\t\t\t\t---------------------------------------\n");
  printf("\t\t\t\t\t     TREE IN PREORDERED ARRAY FORM\n");
  printf("\t\t\t\t\t---------------------------------------\n");
  printf("\n\t\tHIGHEST INTERNAL NODE: %d\t",cst[0][2]);
  printf("LAST TAXON ADDED: %d\t",cst[0][3]);
  printf("COST OF SUB TASK: %d\n\n",cst[0][8]);

  printf("\tLCHILD\tNODE\tRCHILD\tPARENT\t");
  printf("  COST UNTIL HERE\t{STATE VECTOR}\n");

  printf("ROOT NODE\n");
  printf("\t%3d\t root\t%3d\t root\t\t%3d\t\t",cst[0][lhsOffSet+skipSize],     \
		                                    cst[1][1], cst[0][8]);
  printf("{");
  for (ctr=9; ctr < 9+matrix.num_pars_inf_sites; ctr++) {
         printf("%2d ",cst[0][ctr]);
  }
  printf("}\n");

  printf("LEFT TREE\n");
  if (cst[0][lhsOffSet+skipSize] < 0) {
    printf("\t%3d\t%3d\t%3d\t root\t\t%3d\t\t",                                \
		   cst[0][cst[0][lhsOffSet+skipSize+2]],                       \
		   cst[0][lhsOffSet+skipSize],                                 \
		   cst[0][cst[0][lhsOffSet+skipSize+3]],                       \
		   cst[0][lhsOffSet+skipSize+4]);
  } else {
    printf("\t null\t%3d\t null\t root\t\t%3d\t\t", cst[0][lhsOffSet+skipSize],\
		   cst[0][lhsOffSet+skipSize+4]);
  }
  printf("{");
  for (ctr=lhsOffSet+skipSize+5; ctr < lhsOffSet+2*skipSize; ctr++) {
         printf("%2d ",cst[0][ctr]);
  }
  printf("}\n");

  for (ctr=lhsOffSet+2*skipSize; ctr < cst[0][0]; ctr=ctr+skipSize) {
     if (cst[0][ctr] > 0) {
           printf("\t null\t%3d\t null\t%3d\t\t%3d\t\t", cst[0][ctr],          \
			   cst[0][cst[0][ctr+1]], cst[0][ctr+4]);
     } else {
           printf("\t%3d\t%3d\t%3d\t%3d\t\t%3d\t\t",cst[0][cst[0][ctr+2]],     \
	   cst[0][ctr], cst[0][cst[0][ctr+3]], cst[0][cst[0][ctr+1]],          \
	   cst[0][ctr+4]);
     }
     printf("{");
     for (tmp=ctr+5; tmp < ctr+5+matrix.num_pars_inf_sites; tmp++) {
         printf("%2d ",cst[0][tmp]);
     }
     printf("}\n");
  }

  printf("RIGHT TREE\n");
  if (cst[1][1] < 0) {
    printf("\t%3d\t%3d\t%3d\t root\t\t%3d\t\t",cst[1][cst[1][3]],              \
		   cst[1][1], cst[1][cst[1][4]], cst[1][5]);
  } else {
    printf("\t null\t%3d\t null\t root\t\t%3d\t\t", cst[1][1], cst[1][5]);
  }
  printf("{");
  for (ctr=6; ctr < 6+matrix.num_pars_inf_sites; ctr++) {
         printf("%2d ",cst[1][ctr]);
  }
  printf("}\n");

  for (ctr=1+skipSize; ctr < cst[1][0]; ctr=ctr+skipSize) {
     if (cst[1][ctr] > 0) {
           printf("\t null\t%3d\t null\t%3d\t\t%3d\t\t", cst[1][ctr],          \
			   cst[1][cst[1][ctr+1]], cst[1][ctr+4]);
     } else {
           printf("\t%3d\t%3d\t%3d\t%3d\t\t%3d\t\t",cst[1][cst[1][ctr+2]],     \
	   cst[1][ctr], cst[1][cst[1][ctr+3]], cst[1][cst[1][ctr+1]],          \
	   cst[1][ctr+4]);
     }
     printf("{");
     for (tmp=ctr+5; tmp < ctr+5+matrix.num_pars_inf_sites; tmp++) {
         printf("%2d ",cst[1][tmp]);
     }
     printf("}\n");
  }
}
