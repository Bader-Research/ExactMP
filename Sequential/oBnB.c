#include "header.h"

void bnb(void) {
  int level=0, rsize, lsize, cctr, nctr, taxa, cost, indx;

  initBestCost = bestCost;
  numBestSolutions=0;
  //bestCost = 1000053612;

  while (wosac[0][0] != wosac[0][1]) {

     if (level==0) {
       printf("Analyzing task %d\n",wosac[0][0]);
     }
     if (level==matrix.num_taxa-7) {
        for (indx=0; indx < wosac[level][1]; indx++) {
           if (workSpace[level][indx][0][8] < bestCost) {
              bestCost = workSpace[level][indx][0][8];
	      printf("\t\t\t\t\t\t[BEST SCORE: %d]\n",bestCost+pNiCost);
                                                   // Copy LHS Tree
              lsize = workSpace[level][indx][0][0];
              memcpy(solQ[0][0],workSpace[level][indx][0],lsize*sizeof(int)); 
                                                   // Copy RHS tree
              memcpy(solQ[0][1],workSpace[level][indx][1],rsize*sizeof(int));
	      numBestSolutions = 1;
	   } 
#ifdef ALLTREES
	   else if (workSpace[level][indx][0][8] == bestCost) {
                                                   // Copy LHS Tree
                   lsize = workSpace[level][indx][0][0];
                   memcpy(solQ[numBestSolutions][0],workSpace[level][indx][0],lsize*sizeof(int)); 
                                                   // Copy RHS tree
                   rsize = workSpace[level][indx][0][1];
                   memcpy(solQ[numBestSolutions][1],workSpace[level][indx][1],rsize*sizeof(int));
		   ++numBestSolutions;
	   }
#endif
        }
	--level;
     } else {
        cctr = wosac[level][0];
        lsize = workSpace[level][cctr][0][0];
        rsize = workSpace[level][cctr][0][1];
        nctr = 0;
        taxa = taxaQueue[workSpace[level][cctr][0][5]];
        //cost = workSpace[level][cctr][0][8];
        cost = workSpace[level][cctr][0][8] + plbValues[level];
#ifdef ALLTREES
        if (cost <= bestCost) {
#endif

#ifdef ONETREE
        if (cost < bestCost) {
#endif
           for (indx=lhsOffSet+skipSize; indx<=lsize-skipSize; indx+=skipSize) {
              genTaskLeft(workSpace[level][cctr],workSpace[level+1][nctr],taxa,indx);
              nctr = nctr + 1;
           }
           for (indx=1+skipSize; indx <= rsize-skipSize; indx += skipSize) {
              genTaskRight(workSpace[level][cctr],workSpace[level+1][nctr],taxa,indx);
              nctr = nctr + 1;
           }
	   wosac[level+1][0] = 0;
	   ++wosac[level][0];
	   ++level;
	} else {
	   ++wosac[level][0];
	}
     }
     for (cctr=level; cctr >=1; cctr--) {
         if (wosac[cctr][0]==wosac[cctr][1]) {
		 --level;
	 } else {
                 level = cctr;
		 break;
	 }
     }
  }
}
