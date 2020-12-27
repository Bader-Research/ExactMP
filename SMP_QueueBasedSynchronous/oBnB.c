#include "header.h"

void bnb(THREADED) {
  int ****workSpace, *plbValues, ***fitchArray, ***mpTree, **wosac;
  int level=0, rsize, lsize, cctr, nctr, taxa, cost, indx, prevTID=0, currTID=0;
  int done, rem, statesOver, statesRem, site, plbLen, plbIndx, mask,
      bpos, maskRem, maskOver;
  int ctsl, ctsr, ntsl, ntsr, sp=-1, hin=0, rci=0, lci=0, si=0, length=0,
      lc, rc, pa, flag, var, nflg, lim, r, c, ctr, val;
  int lhsSize, rhsSize, row, col, chldCol, pos, mpRows;
  int lhsOffSet, skipSize, *stack;
  int ***localQ, ***tempQ, lqend, tqend, upload, download, baseLevel;
  int *termCond, *copyQ, decompLevel=0, doBnB=0, copied=0, **initTree;
  int *currTask, *prevTask, i, j, sopTask=0, socTask=0, **beginTask, **blankTask;

  double pj;

//on_one_thread {

  skipSize = 5 + matrix.num_pars_inf_sites,
  lhsOffSet = 4;
      
  stack = (int *)calloc(matrix.num_taxa, sizeof(int));
  assert(stack);

  plbValues = (int *)calloc(matrix.num_taxa-6, sizeof(int));
  assert(plbValues);
  for (plbIndx=7; plbIndx < matrix.num_taxa; plbIndx++) {
    for (plbLen=0, site=0; site < matrix.num_pars_inf_sites; site++) {
     for (statesOver=0, done=0; done < plbIndx; done++) {
         statesOver = statesOver |                                             \
                                matrix.reord_sites_enc[taxaQueue[done]-1][site];
     }
     for (statesRem=0, rem=plbIndx; rem < matrix.num_taxa; rem++) {
         statesRem = statesRem |                                               \
                                 matrix.reord_sites_enc[taxaQueue[rem]-1][site];
     }

     mask = 1;
     for (bpos=0; bpos < 32; bpos++) {
         mask = mask << bpos;
         maskRem = statesRem & mask;
         maskOver = statesOver & mask;
         if (maskRem != 0 && maskOver == 0)  {
            ++plbLen;
         }
     }

      plbValues[plbIndx-7] = plbLen;
    }
  }

  fitchArray = (int ***)calloc(16, sizeof(int **));
  assert(fitchArray);
  for (r=0; r < 16; r++) {
     fitchArray[r] = (int **)calloc(16, sizeof(int *));
     assert(fitchArray[r]);
  }
  for (r=0; r < 16; r++) {
     for (c=0; c < 16; c++) {
        fitchArray[r][c] = (int *)calloc(2, sizeof(int));
        assert(fitchArray[r][c]);
     }
  }
  for (r=0; r < 16; r++) {
     for (c=0; c < 16; c++) {
         val = r & c;
         if (val==0) {
                 fitchArray[c][r][0] = r | c;
                 fitchArray[c][r][1] = 1;
         }
         else {
                 fitchArray[c][r][0] = val;
                 fitchArray[c][r][1] = 0;
         }
     }
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

  pj = ((double)jobsPerProc)*((double)THREADS);
  /*
  for (r=8; r < matrix.num_taxa-3; r++) {
    pj = pj/(2*r-5);
    printf("T%d: pj: %f, r: %d\n",MYTHREAD, pj, r);
    fflush(stdout);
    if (pj < 1.0) break;
  }
  decompLevel = r-2;
  */
  for (r=3; r < matrix.num_taxa-3; r++) {
    pj = pj/(2*r-5);
    if (pj < 1.0) break;
  }
  decompLevel = r-2;
  if (decompLevel < 5) decompLevel = 5;

  shQ = (int *)calloc(decompLevel, sizeof(int)); 
  assert(shQ);

  copyQ = (int *)calloc(decompLevel, sizeof(int)); 
  assert(copyQ);

  termCond = (int *)calloc(decompLevel, sizeof(int));
  assert(termCond);
  //decompLevel = 5;
  for (c=0; c < decompLevel; c++) {
   termCond[c] = 2*(c+3)-5;
  }

  beginTask = (int **)calloc(2, sizeof(int *));
  assert(beginTask);
  beginTask[0] = (int *)calloc(2*matrix.num_taxa*skipSize+lhsOffSet, sizeof(int));
  assert(beginTask[0]);
  beginTask[1] = (int *)calloc(2*matrix.num_taxa*skipSize+1, sizeof(int));
  assert(beginTask[1]);

  blankTask = (int **)calloc(2, sizeof(int *));
  assert(blankTask);
  blankTask[0] = (int *)calloc(2*matrix.num_taxa*skipSize+lhsOffSet, sizeof(int));
  assert(blankTask[0]);
  blankTask[1] = (int *)calloc(2*matrix.num_taxa*skipSize+1, sizeof(int));
  assert(blankTask[1]);

  currTask = (int *)calloc(2*matrix.num_taxa, sizeof(int));
  assert(currTask);
  prevTask = (int *)calloc(2*matrix.num_taxa, sizeof(int));
  assert(prevTask);

  /*
  prevTask[0] = 1;
  prevTask[1] = -1;
  prevTask[2] = 2;
  prevTask[3] = 3;
  sopTask = 4;
  decompLevel = 3;
  convert(prevTask, sopTask, beginTask, taxaQueue[decompLevel+1], decompLevel, stack, TH);
  printf("Test Task: ");
  for (i=0; i < socTask; i++) {
    printf(" %d",currTask[i]);
  }
  printf("\n");
  displayTree(beginTask, lhsOffSet, skipSize, TH); 
  getchar();
  */
  on_one_thread {
    printf("Decomposed to level %d\n\n",decompLevel);
    fflush(stdout);
  }
  node_Barrier();

 /* Actual computation starts here */
  while (1) {
   prevTask[0] = 1;
   prevTask[1] = -1;
   prevTask[2] = 2;
   prevTask[3] = 3;
   sopTask = 4;
   if (shQ[0]==1) {
#ifdef debugQ
     printf("T%d: shQ[0]: %d\n",MYHREAD, shQ[0]);
     fflush(stdout);
#endif
     break;
   } else {
      pthread_mutex_lock(&(sqlock));
#ifdef debugQ
      printf("T%d: locking\n", MYTHREAD);
      fflush(stdout);
#endif
      memcpy(copyQ, shQ, decompLevel*sizeof(int));
      memcpy(copyQ, shQ, decompLevel*sizeof(int));
      ++shQ[decompLevel-1];
      for (i=decompLevel-1; i > 0; i--) {
        if (shQ[i]==1+termCond[i]) {
          shQ[i] = 0;
          ++shQ[i-1];
#ifdef debugQ
          printf("T%d: i:%d, incrmenting shQ[%d]: %d\n",MYTHREAD, i, i-1, shQ[i-1]);
#endif
        } else {
          break;
        }
      }
      if (verbose) {
         printf("T%d: working on job # %.0f\n", MYTHREAD, ++totalJobs);
         fflush(stdout);
      }
      pthread_mutex_unlock(&(sqlock));
        doBnB = 0;
        for (r=1; r < decompLevel; r++) {
           taxa = taxaQueue[r+2];
           indx = copyQ[r];
           if (prevTask[indx]!=-1) {
               doBnB = 1;
#ifdef debugQ
               printf("T%d: prevTask[%d] != -1\n",MYTHREAD, indx);
               fflush(stdout);
#endif
               for (j=0; j < indx; j++) {
                 currTask[j] = prevTask[j];
#ifdef debugQ
                 printf("T%d: 1. copying %d of prev into %d of curr\n",MYTHREAD, j,j);
                 fflush(stdout);
#endif
               }
               currTask[indx] = 0-r-1;
               currTask[indx+1] = taxa;
               for (j=indx; j < sopTask; j++) {
#ifdef debugQ
                 printf("T%d: 2. copying %d of prev into %d of curr\n",MYTHREAD, j,j+2);
                 fflush(stdout);
#endif
                 currTask[j+2] = prevTask[j];
               }
               socTask = sopTask+2;
               sopTask = socTask;
               memcpy(prevTask, currTask, socTask*sizeof(int));
#ifdef debugQ
               printf("T%d: 3. memcopied %d of curr into prev\n",MYTHREAD, 2*(r+1)+2);
               fflush(stdout);
#endif
           } else {
               doBnB = 0;
#ifdef debugQ
               printf("T%d: 4. breaking\n", MYTHREAD);
               fflush(stdout);
#endif
               break;
           }
#ifdef debugQ
           if (r==decompLevel-1) {
           printf("Current Task: ");
           for (i=0; i < socTask; i++) {
             printf(" %d",currTask[i]);
           }
           printf("\n");
           }
#endif
        }
        if (doBnB) {
           memcpy(beginTask[0], blankTask[0], (2*matrix.num_taxa*skipSize+lhsOffSet)*sizeof(int));
           memcpy(beginTask[1], blankTask[1], (2*matrix.num_taxa*skipSize+1)*sizeof(int));
           convert(currTask, socTask, beginTask, taxaQueue[decompLevel+1], decompLevel, stack, TH);
#ifdef debugQ
           displayTree(beginTask, lhsOffSet, skipSize, TH); 
#endif

           for (r=0; r <= matrix.num_taxa-7; r++) {
            wosac[r][0] = 0;
           }
           level = abs(beginTask[0][2])-5;            // compute current jobs level
           lsize = beginTask[0][0];
           rsize = beginTask[0][1];
           nctr = 0;
           taxa = taxaQueue[beginTask[0][5]];
           cost = beginTask[0][8] + plbValues[level];

           if (cost <= bestCost) {
             if (level==matrix.num_taxa-7) {
              pthread_mutex_lock(&(sqlock));
              if (cost==bestCost) {
                if (solQCtr < keepTrees) {
                  ++solQCtr;
                  lsize = beginTask[0][0];
                  memcpy(solQ[solQCtr][0],beginTask[0],lsize*sizeof(int));
                                                   /* Copy RHS tree */
                  rsize = beginTask[0][1];
                  memcpy(solQ[solQCtr][1],beginTask[1],rsize*sizeof(int));
                }
              } else if (cost < bestCost) {
                  solQCtr=0;
                  lsize = beginTask[0][0];
                  memcpy(solQ[solQCtr][0],localQ[lqend][0],lsize*sizeof(int));
                                                   /* Copy RHS tree */
                  rsize = localQ[lqend][0][1];
                  memcpy(solQ[solQCtr][1],beginTask[1],rsize*sizeof(int));
                  printf("\t\t\tT%d: Improved best cost to %d\n",MYTHREAD, bestCost+pNiCost);
                  fflush(stdout);
              }
              bestCost = cost;
              pthread_mutex_unlock(&(sqlock));
             } else {
               for (indx=lhsOffSet+skipSize; indx<=lsize-skipSize; indx+=skipSize) {
                 genTaskLeft(beginTask, workSpace[level+1][nctr],taxa,indx,stack,TH);
                 nctr = nctr + 1;
               }
               for (indx=1+skipSize; indx <= rsize-skipSize; indx += skipSize) {
                 genTaskRight(beginTask,workSpace[level+1][nctr],taxa,indx,stack,TH);
                 nctr = nctr + 1;
               }
               loadDistribution[MYTHREAD] = loadDistribution[MYTHREAD] + (double)(nctr);
               level = level+1;
               baseLevel = level;
               for (ctr=1; ctr <= matrix.num_taxa-7; ctr++) {
                wosac[ctr][0]=0;
               }
               currTID = 0;
               prevTID = 0;
               while (wosac[baseLevel][0] <= wosac[baseLevel][1]) {
                prevTID = currTID;
                currTID = wosac[baseLevel][0];
                if (prevTID==wosac[baseLevel][1] && currTID==0) break;
#ifdef debug
                   printf("\nStatus at last level: ");
                   fflush(stdout);
                   for (ctr=0; ctr < matrix.num_taxa-7; ctr++) {
                      printf("wosac[%d][0]: %d   ",ctr, wosac[ctr][0]);
                      fflush(stdout);
                   }
#endif
                   if (level==matrix.num_taxa-7) {
#ifdef debug
                     printf("\nwosac[%d][0]:",level);
                     fflush(stdout);
#endif
                     for (wosac[level][0]=0; wosac[level][0] < wosac[level][1]; wosac[level][0]++) {
#ifdef debug
                      printf("  %d",wosac[level][0]);
                      fflush(stdout);
#endif
                      if (workSpace[level][wosac[level][0]][0][8] < bestCost) {
                         pthread_mutex_lock(&(bclock));
                         bestCost = workSpace[level][wosac[level][0]][0][8];
	                 printf("\t\t\t\t\t\tT%d: [BEST SCORE: %d]\n",MYTHREAD, bestCost+pNiCost);
                         fflush(stdout);
                                                   // Copy Tree to Solution Queue
                         lsize = workSpace[level][wosac[level][0]][0][0];
                         rsize = workSpace[level][wosac[level][0]][0][1];
                         memcpy(solQ[0][0],workSpace[level][wosac[level][0]][0],lsize*sizeof(int));
                         memcpy(solQ[0][1],workSpace[level][wosac[level][0]][1],rsize*sizeof(int));
	                 solQCtr=0;
                         pthread_mutex_unlock(&(bclock));
	              }
                      else if (workSpace[level][wosac[level][0]][0][8] == bestCost) {
                         pthread_mutex_lock(&(bclock));
                         if (solQCtr < keepTrees-1) {
                                                   // Copy Tree to Solution Queue
	                   ++solQCtr;
                           lsize = workSpace[level][wosac[level][0]][0][0];
                           rsize = workSpace[level][wosac[level][0]][0][1];
                           memcpy(solQ[solQCtr][0],workSpace[level][wosac[level][0]][0],lsize*sizeof(int));
                           memcpy(solQ[solQCtr][1],workSpace[level][wosac[level][0]][1],rsize*sizeof(int));
                         }
                         pthread_mutex_unlock(&(bclock));
	              }
                     }
                     loadDistribution[MYTHREAD] += wosac[level][1];
#ifdef debug
                     printf("\n");
                     fflush(stdout);
#endif
                } else {
                   nctr = 0;
                   cctr = wosac[level][0];
                   lsize = workSpace[level][cctr][0][0];
                   rsize = workSpace[level][cctr][0][1];
                   taxa = taxaQueue[workSpace[level][cctr][0][5]];
                   cost = workSpace[level][cctr][0][8] + plbValues[level];
                   if (cost <= bestCost) {
                    for (indx=lhsOffSet+skipSize; indx<=lsize-skipSize; indx+=skipSize) {
                       genTaskLeft(workSpace[level][cctr],workSpace[level+1][nctr],taxa,indx,stack,TH);
                       nctr = nctr + 1;
                    }
                    for (indx=1+skipSize; indx <= rsize-skipSize; indx += skipSize) {
                       genTaskRight(workSpace[level][cctr],workSpace[level+1][nctr],taxa,indx,stack,TH);
                       nctr = nctr + 1;
                    }
	            wosac[level+1][0] = 0;
	            ++wosac[level][0];
	            ++level;
                    loadDistribution[MYTHREAD] += nctr;
	           } else {
	            ++wosac[level][0];
                    loadDistribution[MYTHREAD] += 1;
	           }
                }
                for (cctr=level; cctr >=baseLevel; cctr--) {
                  if (wosac[cctr][0]==wosac[cctr][1]) {
		     --level;
	          } else {
                     level = cctr;
		     break;
	          }
                }
                if (level==baseLevel-1) break;
#ifdef debug
                    printf("\n1-l: Wrapping up to level %d",level);
                    fflush(stdout);
#endif
                }
             }
           }
        }
   }
  }
//}
}
