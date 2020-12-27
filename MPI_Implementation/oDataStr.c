#include "header.h"

void allocFitchArray(void) {

  int lc, rc, val;

  fitchArray = (int ***)calloc(16, sizeof(int **));
  assert(fitchArray);
  for (rc=0; rc < 16; rc++) {
     fitchArray[rc] = (int **)calloc(16, sizeof(int *));
     assert(fitchArray[rc]);
  }
  for (rc=0; rc < 16; rc++) {
     for (lc=0; lc < 16; lc++) {
        fitchArray[rc][lc] = (int *)calloc(2, sizeof(int));
        assert(fitchArray[rc][lc]);
     }
  }

  for (rc=0; rc < 16; rc++) {
     for (lc=0; lc < 16; lc++) {
         val = rc & lc;
         if (val==0) {
		 fitchArray[lc][rc][0] = rc | lc;
		 fitchArray[lc][rc][1] = 1;
	 }
	 else {
		 fitchArray[lc][rc][0] = val;
		 fitchArray[lc][rc][1] = 0;
	 }
     }
  }

  /*
  printf("\n");
  for (rc=1; rc < 16; rc++) {
     for (lc=1; lc < 16; lc++) {
         printf("rc: %d, lc: %d, pa: %d, len: %d\n", rc, lc, fitchArray[lc][rc][0], fitchArray[lc][rc][1]);
     }
     printf("\n");
  }
  */

}

void allocMem4WorkSpace(void) {
  int r, c, lim, add;

                                        // allocate memory for the work space
					// this includes all the levels of the 
					// mp tree starting from row with 7 
                                        // taxa. However, other than row 0 in 
                                        // the work space, each row contains
					// only (2*level-5) trees can be 
					// allocated (done to keep memory
					// under control). Row 0 contains all
                                        // (2*7-5)!! = 945 trees
  lim = matrix.num_taxa-6;
  workSpace = (int ****)calloc(lim, sizeof(int ***));
  assert(workSpace);
  workSpace[0] = (int ***)calloc(945, sizeof(int **));
  assert(workSpace[0]);
  for (c=0; c < 945; c++) {
      workSpace[0][c] = (int **)calloc(2, sizeof(int *));
      workSpace[0][c][0] = (int *)calloc(12*skipSize+lhsOffSet, sizeof(int));
      workSpace[0][c][1] = (int *)calloc(12*skipSize+1, sizeof(int));
  }

  for (r=1; r < lim; r++) {
      workSpace[r] = (int ***)calloc((2*(r+7)-5), sizeof(int **));
      for (c=0; c < 2*(r+7)-5; c++) {
          workSpace[r][c] = (int **)calloc(2, sizeof(int *));
          workSpace[r][c][0] = (int *)calloc((2*(r+7)-1)*skipSize+lhsOffSet,   \
			                                           sizeof(int));
          workSpace[r][c][1] = (int *)calloc((2*(r+7)-1)*skipSize+1,           \
			                                           sizeof(int));
      }
  }

  solQ = (int ***)calloc(keepTrees, sizeof(int **));
  for (c=0; c < keepTrees; c++) {
      solQ[c] = (int **)calloc(2, sizeof(int *));
      solQ[c][0] = (int *)calloc((2*matrix.num_taxa-1)*skipSize+lhsOffSet,     \
		                                                   sizeof(int));
      solQ[c][1] = (int *)calloc((2*matrix.num_taxa-1)*skipSize+1, sizeof(int));
  }

  lim = matrix.num_taxa-6;
  wosac = (int **)calloc(lim, sizeof(int *));
  for (r=0; r < lim; r++) {
      wosac[r] = (int *)calloc(2, sizeof(int));
  }

  wosac[0][1] = 945;
  for (r=1; r < lim; r++) {
      wosac[r][1] = 2*(r+7)-5;
  }

  inlev = (int *)calloc(numProcs-1, sizeof(int));
  mywoli = (int **)calloc(numProcs-1, sizeof(int *));
  for (r=0; r < numProcs-1; r++) {
    mywoli[r] = (int *)calloc(2,sizeof(int));
    inlev[r] = -1;
  }

  if (numProcs <= 945) {
     lim = 945/(numProcs-1);
     mywoli[0][0] = 0;
     mywoli[0][1] = lim;
     inlev[0] = 0;
     for (r=1; r < numProcs-1; r++) {
       inlev[r] = 0;
       mywoli[r][0] = mywoli[r-1][1];
       mywoli[r][1] = mywoli[r][0]+lim;
     }
  }

  lim = 945 % (numProcs-1);
  if (lim > 0) {
     add=1;
     for (r=numProcs-lim-1; r < numProcs-1; r++) {
       printf("r: %d, add: %d\n",r,add);
       fflush(stdout);
       mywoli[r][0] = mywoli[r-1][1];
       mywoli[r][1] = mywoli[r][1]+add;
       ++add;
     }
  }

  if (my_rank==0) {
     printf("\nWork Space Distribution: \n");
     fflush(stdout);
     for (r=0; r < numProcs-1; r++) {
        printf("\tRank %3d: \t%3d\t-\t%3d\n",r,mywoli[r][0],mywoli[r][1]);
        fflush(stdout);
     }
     for (r=0; r < numProcs-1; r++) {
        printf("\tInitial level for %3d: %3d\n",r,inlev[r]);
        fflush(stdout);
     }
  }

  if (my_rank!=numProcs-1) {
    my_level = inlev[my_rank];
    wosac[0][0] = mywoli[my_rank][0];
    wosac[0][1] = mywoli[my_rank][1];
  }
}

void getKthLevel(int mpRows) {
 int ***mpTree, lhsSize, rhsSize,
     row, col, chldCol, pos, ctr, sp,
     rci=0, lci=0, si=0, length=0, site=0, lc, rc, pa, flag, var, nflg;

 mpTree = (int ***)calloc(mpRows+1, sizeof(int **));
 
                                        // This piece of code just generates
					// all trees at kth level, where 
					// k=mpRows. This algorithm is very 
					// inefficient though, but it is used 
					// only once in the program...
					// It recursively finds all trees in
					// next level to the current level
					// until the requested level is reached
 for (row=1; row <= mpRows; row++) {
   mpTree[row] = (int **)calloc(numNodesInRow(row), sizeof(int *));
   for (col=0; col < numNodesInRow(row); col++) {
       mpTree[row][col] = (int *)calloc(2*(row+1), sizeof(int ));
   }
 }

 mpTree[1][0] = (int *)calloc(4, sizeof(int ));
 mpTree[1][0][0] =  1;
 mpTree[1][0][1] = -1;
 mpTree[1][0][2] =  2;
 mpTree[1][0][3] =  3;

 for (row=1; row < mpRows; row++) {
   for (chldCol=0, col=0; col < numNodesInRow(row); col++) {
      for (pos=0; pos < 2*(row+1); pos++) {
         if (mpTree[row][col][pos] != -1) {
            for (ctr=0; ctr < pos; ctr++) {
               mpTree[row+1][chldCol][ctr] = mpTree[row][col][ctr];
	    }
            mpTree[row+1][chldCol][pos] = -row-1;
            mpTree[row+1][chldCol][pos+1] = taxaQueue[row+2];
            for (ctr=pos; ctr < 2*(row+1); ctr++) {
               mpTree[row+1][chldCol][ctr+2] = mpTree[row][col][ctr];
	    }
	    chldCol += 1;
	 }
      }
   } 
 }
					// Here, all the trees until requested
					// level are printed out
 /*
 printf("\n");
 for (row=1; row <= mpRows; row++) {
   chldCol = 0;
   for (col=0; col < numNodesInRow(row); col++) {
	   if (row!=1 && col%(2*row-1) == 0) {
		   chldCol += 3;
	   }
	   for (ctr=0; ctr < chldCol; ctr++) {
		    printf(" ");
	   }
	   printf("[%d][%d]: ",row,col);
      for (pos=0; pos < 2*(row+1); pos++) {
            printf(" %d", mpTree[row][col][pos]);
      }
      printf("\n");
   }
      printf("\n\n");
 } */

					// For all nodes in the current level
					// do the following-
					// (a) set important fields (in offset)
					// (b) copy states of leaf nodes
					// (c) set parent-child relationships
					// (d) compute cost of each of 945 trees
 for (col=0; col < numNodesInRow(mpRows); col++) {
     lhsSize = 0;
     for (pos = 0; pos < 2*(mpRows+1); pos++) {
         if (mpTree[5][col][pos] != -1) {
            ++lhsSize;
         } else {
            chldCol = pos;
            break;
	 }
     }
     rhsSize = 0;
     for (pos = chldCol; pos < 2*(mpRows+1); pos++) {
         ++rhsSize;
     }
                                        // SET IMPORTANT FIELDS OF LHS TREE
                                        // Size of LHS tree
     workSpace[0][col][0][0] = (lhsSize+1)*skipSize + lhsOffSet;  
                                        // Size of RHS tree
     workSpace[0][col][0][1] = rhsSize*skipSize + 1;       
                                        // Highest Internal Node
     workSpace[0][col][0][2] = -5;
                                        // Last taxa added
     workSpace[0][col][0][3] = taxaQueue[5];
                                        // dummy
     workSpace[0][col][0][4] = ROOT;                       
                                        // Next position in taxa Queue
     workSpace[0][col][0][5] = 7;
                                        // cost of task
     workSpace[0][col][0][8] = 0;                          
                                        // parent of lhs root is ROOT
     workSpace[0][col][0][lhsOffSet+skipSize+1] = 4;          
                                        // parent of rhs root is ROOT
     workSpace[0][col][1][2] = 4;                          

                                        // SET IMPORTANT FIELDS OF RHS TREE
                                        // Size of RHS tree
     workSpace[0][col][1][0] = rhsSize*skipSize + 1;       

                                        // Convert 1D array into 2D LHS/RHS
					// array with states loaded
					// This for loop covers LHS tree and 
					// breaks at RHS tree
     for (pos = 0; pos < 2*(mpRows+1); pos++) {
                                        // Get position where RHS tree starts
         if (mpTree[5][col][pos] == -1) {
		 chldCol = pos;
		 break;
	 }

                                        // copy taxon ID
	 workSpace[0][col][0][(pos+1)*skipSize+lhsOffSet] = mpTree[5][col][pos];
                                        // if leaf node, copy state vector
	 if (mpTree[5][col][pos] > 0) {
            memcpy(&workSpace[0][col][0][(pos+1)*skipSize+lhsOffSet+5],        \
			    matrix.reord_sites_enc[mpTree[5][col][pos]-1],     \
			    matrix.num_pars_inf_sites*sizeof(int));
	    workSpace[0][col][0][((pos+1)*skipSize)+lhsOffSet+4] = 0;
	 }
     }

					// This for loop covers RHS tree
     for (row = chldCol; row < 2*(mpRows+1); row++) {
                                        // copy taxon ID
	 workSpace[0][col][1][(row-chldCol)*skipSize+1] = mpTree[5][col][row];
                                        // if leaf node, copy state vector
	 if (mpTree[5][col][row] > 0) {
            memcpy(&workSpace[0][col][1][(row-chldCol)*skipSize+1+5],          \
			    matrix.reord_sites_enc[mpTree[5][col][row]-1],     \
			    matrix.num_pars_inf_sites*sizeof(int));
	    workSpace[0][col][1][(row-chldCol)*skipSize+1+4] = 0;
	 }
     }
 }

                                        // set LHS tree parent child relation
 for (col=0; col < numNodesInRow(mpRows); col++) {
      row = workSpace[0][col][0][0];
      for (sp=-1,ctr=row-skipSize;ctr>=lhsOffSet+skipSize;ctr=ctr-skipSize) {
       if (workSpace[0][col][0][ctr] > 0) {
         stack[++sp] = ctr;
       } else {
         workSpace[0][col][0][stack[sp]+1] = ctr;
	 workSpace[0][col][0][ctr+2] = stack[sp];

	 sp = sp-1;
	 workSpace[0][col][0][stack[sp]+1] = ctr;
	 workSpace[0][col][0][ctr+3] = stack[sp];

	 stack[sp] = ctr;
       }
      }

                                        // set RHS tree parent child relation
      row = workSpace[0][col][1][0];
      for (sp=-1, ctr=row-skipSize; ctr > 0; ctr = ctr-skipSize) {
       if (workSpace[0][col][1][ctr] > 0) {
         stack[++sp] = ctr;
       } else { 
         workSpace[0][col][1][stack[sp]+1] = ctr;
	 workSpace[0][col][1][ctr+2] = stack[sp];

	 sp = sp-1;
	 workSpace[0][col][1][stack[sp]+1] = ctr;
	 workSpace[0][col][1][ctr+3] = stack[sp];

	 stack[sp] = ctr;
       }
      }
      
                                        // compute LHS cost
  for (si=workSpace[0][col][0][0]-skipSize; si > 4; si = si - skipSize) {
        if (workSpace[0][col][0][si] < 0) {
          lci = workSpace[0][col][0][si+2];
          rci = workSpace[0][col][0][si+3];
          for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {
             lc = workSpace[0][col][0][lci+site+5];
	     rc = workSpace[0][col][0][rci+5+site];
             pa = lc | rc;
             var = lc & rc;
             flag = !var;
             nflg = ~flag;
             nflg = nflg + 1;
             pa = pa & nflg;
             pa = pa | var;
             length = length+flag;
             workSpace[0][col][0][si+5+site] = pa;
          }
          workSpace[0][col][0][si+4] = workSpace[0][col][0][lci+4]+workSpace[0][col][0][rci+4]+length;
        }
  }

                                        // compute RHS cost
  for (si=workSpace[0][col][1][0]-skipSize; si >= 1; si = si - skipSize) {
        if (workSpace[0][col][1][si] < 0) {
          lci = workSpace[0][col][1][si+2];
          rci = workSpace[0][col][1][si+3];
          for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {
             lc = workSpace[0][col][1][lci+site+5];
	     rc = workSpace[0][col][1][rci+5+site];
             pa = lc | rc;
             var = lc & rc;
             flag = !var;
             nflg = ~flag;
             nflg = nflg + 1;
             pa = pa & nflg;
             pa = pa | var;
             length = length+flag;
             workSpace[0][col][1][si+5+site] = pa;
          }
          workSpace[0][col][1][si+4] = workSpace[0][col][1][lci+4]+workSpace[0][col][1][rci+4]+length;
        }
  }

                                        // compute ROOT cost
  for (length=0, site=0; site < matrix.num_pars_inf_sites; site++) {
     lc = workSpace[0][col][0][lhsOffSet+skipSize+5+site];
     rc = workSpace[0][col][1][6+site];
     pa = lc | rc;
     var = lc & rc;
     flag = !var;
     nflg = ~flag;
     nflg = nflg + 1;
     pa = pa & nflg;
     pa = pa | var;
     length = length+flag;
     workSpace[0][col][0][9+site] = pa;
  }
  workSpace[0][col][0][8] = workSpace[0][col][0][lhsOffSet+skipSize+4]+workSpace[0][col][1][5]+length;
 }

 for (row=1; row < mpRows; row++) {
   for (col=0; col < numNodesInRow(row); col++) {
       free(mpTree[row][col]);
   }
   free(mpTree[row]);
 }
 
 for (col=0; col < numNodesInRow(mpRows); col++) {
       free(mpTree[mpRows][col]);
 }
 free(mpTree[row]);
 free(mpTree);

}

int numNodesInRow(int row) {
  if (row == 1) return 1;
  else return (2*row-1)*numNodesInRow(row-1);
}

void clearWorkSpace(void) {
  int r, c, lim;

  for (c=0; c < 945; c++) {
      free(workSpace[0][c][0]);
      free(workSpace[0][c][1]);
      free(workSpace[0][c]);
  }
  free(workSpace[0]);

  lim = matrix.num_taxa-6;
  for (r=1; r < lim; r++) {
      for (c=0; c < 2*(r+7)-5; c++) {
          free(workSpace[r][c][0]);
          free(workSpace[r][c][1]);
          free(workSpace[r][c]);
      }
      free(workSpace[r]);
  }
  free(workSpace);

  for (c=0; c < keepTrees; c++) {
      free(solQ[c][0]);
      free(solQ[c][1]);
      free(solQ[c]);
  }
  free(solQ);

  lim = matrix.num_taxa-6;
  for (r=0; r < lim; r++) {
      free(wosac[r]);
  }
  free(wosac);

  free(stack);
  free(taxaQueue);
  free(plbValues);

  for (r=0; r < 16; r++) {
     for (c=0; c < 16; c++) {
        free(fitchArray[r][c]);
     }
  }
  for (r=0; r < 16; r++) {
     free(fitchArray[r]);
  }
  free(fitchArray);

}
