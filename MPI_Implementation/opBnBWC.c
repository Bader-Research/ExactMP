#include "header.h"

void paBnB(void) {
  int level=0, rsize, lsize, cctr, nctr, taxa, cost, indx, mnr=0;

                                           // If I am the last processor, I
                                           // wouldnt do crap work :)
  if (my_rank==numProcs-1) {
       while (1) {
         flag = 0;
         MPI_Iprobe(MPI_ANY_SOURCE,                                           \ 
                                 2,                                           \
                    MPI_COMM_WORLD,                                           \
                             &flag,                                           \
                          &status);
         if (flag==1) {
            MPI_Irecv(&nBC,                                                   \
                         1,                                                   \
                           MPI_INT,                                           \
                    MPI_ANY_SOURCE,                                           \
                                 2,                                           \
                    MPI_COMM_WORLD,                                           \
                         &request);
            MPI_Wait(&request, &status);
            ++mnr;
         }
         if (mnr==numProcs-1) {
            for (cctr=0; cctr < numProcs-1; cctr++) {
               MPI_Isend(&my_rank,                                           \
                                1,                                           \
                          MPI_INT,                                           \
                             cctr,                                           \
                                1,                                           \
                   MPI_COMM_WORLD,                                           \
                        &request);
               MPI_Wait(&request, &status);
            }
            goto TERMINATE;
         }
       }
  } else {
                                           // Else, I do the main job of   
                                           // computing the global minima
       //printf("rank: %d, else block step 1\n",my_rank);
       //fflush(stdout);
       initBestCost = bestCost;
       numBestSolutions=0;
       //bestCost = 600;

                                           // Only if I have a job at hand,
                                           // will I get in...
       if (my_level >= 0) {
LOOP1: 
          //printf("rank: %d, else block LOOP1 step 2\n",my_rank);
          //fflush(stdout);
          level = my_level;
          mnr = 0;                         // my negative replies = 0 as of now

                                           // as long as I have jobs at hand, I
                                           // slog on...
          while (wosac[my_level][0] < wosac[my_level][1]) {
                                           // record current best cost to 
                                           // propagate if it improves
              oBC = bestCost;
              //printf("rank: %d, else block LOOP1 step 3, wosac[%d][0]: %d, wosac[%d][1]: %d, level: %d\n",my_rank,my_level,wosac[my_level][0], my_level, wosac[my_level][1], level);
              //fflush(stdout);
                                           // if this is the last level,  
                                           // check for best cost change
              if (level==matrix.num_taxa-7) {
                //printf("rank: %d, else block LOOP1 step 4\n",my_rank);
                //fflush(stdout);
                for (indx=0; indx < wosac[level][1]; indx++) {
                 if (workSpace[level][indx][0][8] < bestCost) {
                    //printf("rank: %d, else block LOOP1 step 5\n",my_rank);
                    //fflush(stdout);
                    bestCost = workSpace[level][indx][0][8];
	            printf("\t\t\t\t\t\tRank[%d]: [BEST SCORE: %d]\n",my_rank, bestCost+pNiCost);
                    fflush(stdout);
                                                   // Copy LHS Tree
                    //lsize = workSpace[level][indx][0][0];
                    //memcpy(solQ[0][0], workSpace[level][indx][0],             \
                                       lsize*sizeof(int)); 
                                                   // Copy RHS tree
                    //rsize = workSpace[level][indx][0][1];
                    //memcpy(solQ[0][1], workSpace[level][indx][1],             \
                                       rsize*sizeof(int));
      	            numBestSolutions = 1;
	         } 
                }
        	wosac[level][0]=wosac[level][1];
                if (level > my_level) {
                     //printf("rank: %d, my level: %d, decrementing level from %d to %d LOOP1 step 44\n",my_rank, my_level, level, level-1);
                     //fflush(stdout);
	             --level;
                }
              } else {
                                           // if this is not last level,  
                                           // expand job into next level & loop
                //printf("rank: %d, else block LOOP1 step 6\n",my_rank);
                //fflush(stdout);
                cctr = wosac[level][0];
                   //printf("rank: %d, else block LOOP1 step 6-A cctr:%d, level: %d, my level: %d\n",my_rank,cctr, level, my_level);
                   //fflush(stdout);
                lsize = workSpace[level][cctr][0][0];
                   //printf("rank: %d, else block LOOP1 step 6-B\n",my_rank);
                   //fflush(stdout);
                rsize = workSpace[level][cctr][0][1];
                   //printf("rank: %d, else block LOOP1 step 6-C\n",my_rank);
                   //fflush(stdout);
                nctr = 0;
                taxa = taxaQueue[workSpace[level][cctr][0][5]];
                   //printf("rank: %d, else block LOOP1 step 6-D\n",my_rank);
                   //fflush(stdout);
                //cost = workSpace[level][cctr][0][8];
                cost = workSpace[level][cctr][0][8] + plbValues[level];
                   //printf("rank: %d, else block LOOP1 step 6-E\n",my_rank);
                   //fflush(stdout);
                if (cost < bestCost) {
                   //printf("rank: %d, else block LOOP1 step 7\n",my_rank);
                   //fflush(stdout);
                   for (indx=lhsOffSet+skipSize; indx<=lsize-skipSize;        \
                                                            indx+=skipSize) {
                      genTaskLeft(workSpace[level][cctr],                     \
                                      workSpace[level+1][nctr],taxa,indx);
                      nctr = nctr + 1;
                   }
                   for (indx=1+skipSize; indx <= rsize-skipSize;              \
                                                          indx += skipSize) {
                      genTaskRight(workSpace[level][cctr],                    \
                                         workSpace[level+1][nctr],taxa,indx);
                      nctr = nctr + 1;
                   }
        	   wosac[level+1][0] = 0;
        	   ++wosac[level][0];
        	   ++level;
        	} else {
                   //printf("rank: %d, else block LOOP1 step 8\n",my_rank);
                   //fflush(stdout);
        	   ++wosac[level][0];
        	}
               }
               //printf("rank: %d, else block LOOP1 step 9\n",my_rank);
               //fflush(stdout);
                                           // check if this job was the last job
               for (cctr=level; cctr > my_level; cctr--) {
                 if (wosac[level][0]==wosac[level][1]) {
                         //printf("rank: %d, my level: %d, decrementing level from %d to %d LOOP1 step 55\n",my_rank, my_level, level, level-1);
                         //fflush(stdout);
	        	 --level;
        	 } else {
        		 break;
        	 }
               }
               //printf("rank: %d, else block LOOP1 step 10\n",my_rank);
               //fflush(stdout);

                    // PARALLEL STUFF: check messages for me
                                           // if better best cost, propagate it
                                           // to all buddies, but not to self
               if (bestCost < oBC) {
                   for (dest=0; dest < numProcs-1; dest++) {
                      if (dest!=my_rank) {
			 //printf("%d sending best cost of %d to %d\n",my_rank, bestCost, dest);
                         //fflush(stdout);
                         MPI_Isend(&bestCost,1,MPI_INT,dest,                  \
                                      3,MPI_COMM_WORLD,&request);
                         MPI_Wait(&request, &status);
                      }
                   }
               }
               //printf("rank: %d, else block LOOP1 step 11\n",my_rank);
               //fflush(stdout);
                                           // if any incoming messages are here
                                           // finish them, before proceeding
               while (1) {
                                           // tag=3 : better-bestcost message
                   flag = 0;
                   MPI_Iprobe(MPI_ANY_SOURCE,3,MPI_COMM_WORLD,&flag,&status);
                   if (flag==1) {
                      MPI_Irecv(&nBC, 1, MPI_INT, MPI_ANY_SOURCE,             \
                                3, MPI_COMM_WORLD, &request);
                      MPI_Wait(&request, &status);
		      //printf("%d receiving best cost of %d from %d\n",my_rank, nBC, status.MPI_SOURCE);
                      //fflush(stdout);
                      if (nBC < bestCost) {
			 //printf("%d updating best cost of %d with %d\n",my_rank, bestCost, nBC);
                         //fflush(stdout);
                         bestCost = nBC;
                      }
                   }
                   //printf("rank: %d, else block LOOP1 step 12-A\n",my_rank);
                   //fflush(stdout);
                                           // tag=4 : job request message
                   flag = 0;
                   MPI_Iprobe(MPI_ANY_SOURCE,4,MPI_COMM_WORLD,&flag,&status);
                   //printf("rank: %d, else block LOOP1 step 12-B\n",my_rank);
                   //fflush(stdout);
                   if (flag==1) {
                      MPI_Irecv(&dest, 1, MPI_INT, MPI_ANY_SOURCE,
                                4, MPI_COMM_WORLD, &request);
                      MPI_Wait(&request, &status);
                      //printf("rank: %d, else block LOOP1 step 13\n",my_rank);
                      //fflush(stdout);
                                           // if jobs at hand, pack and send
                      if (wosac[my_level][0] < wosac[my_level][1] &&          \
                                                my_level < matrix.num_taxa-7) {
                         //printf("rank: %d, else block LOOP1 step 14-A\n",my_rank);
                         //fflush(stdout);
                         indx = wosac[my_level][0];
                         //printf("rank: %d, else block LOOP1 step 14-B\n",my_rank);
                         //fflush(stdout);
                         lsize = workSpace[my_level][indx][0][0];
                         //printf("rank: %d, else block LOOP1 step 14-C\n",my_rank);
                         //fflush(stdout);
                         rsize = workSpace[my_level][indx][0][1];
                         //printf("rank: %d, else block LOOP1 step 14-D\n",my_rank);
                         //fflush(stdout);
                         job[0] = my_level;
                         //printf("rank: %d, else block LOOP1 step 14-E\n",my_rank);
                         //fflush(stdout);
                         memcpy(&job[1], workSpace[my_level][indx][0],        \
                                                         lsize*sizeof(int));
                         //printf("rank: %d, else block LOOP1 step 14-F\n",my_rank);
                         //fflush(stdout);
                         memcpy(&job[1+lsize], workSpace[my_level][indx][1],  \
                                                         rsize*sizeof(int));
                         //printf("rank: %d, else block LOOP1 step 14-G\n",my_rank);
                         //fflush(stdout);
                                           // tag=5 : job available and sending
                         MPI_Isend(job, sopam, MPI_INT, dest,                 \
                                   5, MPI_COMM_WORLD, &request);
                         MPI_Wait(&request, &status);
                         ++wosac[my_level][0];
                         //printf("rank: %d, else block LOOP1 step 15 sent a job to %d\n",my_rank,dest);
                         //fflush(stdout);
                      }
                      else {
                                           // if no jobs at hand, send neg reply
                                           // tag=6
                         //printf("rank: %d, else block LOOP1 step 16\n",my_rank);
                         //fflush(stdout);
                         MPI_Isend(&my_rank, 1, MPI_INT, dest,                \
                                   6, MPI_COMM_WORLD, &request);
                         MPI_Wait(&request, &status);
                      }
                         //printf("rank: %d, else block LOOP1 step 17\n",my_rank);
                         //fflush(stdout);
                   }
                   else if (flag==0) {
                                           // if no messages, move on...
                      //printf("rank: %d, else block LOOP1 step 18\n",my_rank);
                      //fflush(stdout);
                      break;
                   }
               } // WHILE (1)
             } // UNTIL HERE IT WAS BNB PROCESSING
           } // if (my_level >= 0)
           //printf("rank: %d, else block LOOP1 step 14\n",my_rank);
           //fflush(stdout);

                    // PARALLEL STUFF: loop here after retirement
        while (1) {
LOOP2:
                                           // probe & act on terminate message
           flag = 0;
           MPI_Iprobe(numProcs-1,1,MPI_COMM_WORLD,&flag,&status);
           if (flag==1) {
              //printf("rank: %d, LOOP2 step 1-A\n",my_rank);
              //fflush(stdout);
              MPI_Irecv(&dest, 1, MPI_INT, numProcs-1,                         \
                        1, MPI_COMM_WORLD, &request);
              MPI_Wait(&request, &status);
              goto TERMINATE;
           }
           //printf("rank: %d, LOOP2 step 1-B\n",my_rank);
           //fflush(stdout);
                                           // probe & reply job request message
           flag = 0;
           MPI_Iprobe(MPI_ANY_SOURCE,4,MPI_COMM_WORLD,&flag,&status);
           if (flag==1) {
              //printf("rank: %d, LOOP2 step 2-A\n",my_rank);
              //fflush(stdout);
              MPI_Irecv(&dest, 1, MPI_INT, MPI_ANY_SOURCE,                    \
                        4, MPI_COMM_WORLD, &request);
              MPI_Wait(&request, &status);
                           // SEND NO JOB MESSAGE
              MPI_Isend(&my_rank, 1, MPI_INT, dest,                           \
                        6, MPI_COMM_WORLD, &request);
              MPI_Wait(&request, &status);
           }
           //printf("rank: %d, LOOP2 step 2-B\n",my_rank);
           //fflush(stdout);
                           // PROBE FOR BEST COST MESSAGE
                                           // probe & act bestcost message
           flag = 0;
           MPI_Iprobe(MPI_ANY_SOURCE,3,MPI_COMM_WORLD,&flag,&status);
           if (flag==1) {
              //printf("rank: %d, LOOP2 step 3-A\n",my_rank);
              //fflush(stdout);
              MPI_Irecv(&nBC, 1, MPI_INT, MPI_ANY_SOURCE,                     \
                        3, MPI_COMM_WORLD, &request);
              MPI_Wait(&request, &status);
	      //printf("%d receiving best cost of %d from %d\n",my_rank, nBC, status.MPI_SOURCE);
              //fflush(stdout);
              if (nBC < bestCost) {
                         bestCost = nBC;
                         //printf("%d updating best cost of %d with %d\n",my_rank, bestCost, nBC);
                         //fflush(stdout);
              }
           }
           //printf("rank: %d, LOOP2 step 3-B\n",my_rank);
           //fflush(stdout);

                                           // if all buddies done, wait on
                                           // terminate message from slughead
           if (mnr==numProcs-1) {
              //printf("rank: %d, LOOP2 step 4-A\n",my_rank);
              //fflush(stdout);
              goto LOOP2; 
           }
                                           // if no messages, and more jobs
                                           // available, move on...
           if (flag==0) {
              //printf("rank: %d, LOOP2 step 4-B\n",my_rank);
              //fflush(stdout);
               break;
           }
        }
        //printf("rank: %d, LOOP2 step 4-C\n",my_rank);
        //fflush(stdout);

                    // PARALLEL STUFF: request a job, if none wait to terminate
LOOP3:
                                           // goto neighbor first, loop if reqd
           if ((dest=my_rank+1)==numProcs-1) dest=0;
           //printf("rank: %d, LOOP3 step 1\n",my_rank);
           //fflush(stdout);
                                           // request all buddies for job
           for (mnr=1, nbr=0; nbr < numProcs-1; nbr++) {
             //printf("rank: %d, LOOP3 step 1-0\n",my_rank);
             //fflush(stdout);
             if (nbr!=my_rank) {
                flag = 0;
                                           // first check for any incoming 
                                           // job requests & send neg replies
                MPI_Iprobe(dest,4,MPI_COMM_WORLD,&flag,&status);
                if (flag==1) {
                   MPI_Irecv(&nBC, 1, MPI_INT, dest,                          \
                              4, MPI_COMM_WORLD, &request);
                   MPI_Wait(&request, &status);
                   MPI_Isend(&my_rank,1,MPI_INT,dest,6,MPI_COMM_WORLD,&request);
                   MPI_Wait(&request,&status);
                   ++mnr;
                                           // if no job in system, wait to 
                                           // terminate
                   if (mnr==numProcs-1) {
                       //printf("rank: %d, LOOP3 step 1-2\n",my_rank);
                       //fflush(stdout);
                       MPI_Isend(&my_rank, 1, MPI_INT, numProcs-1,            \
                                 2, MPI_COMM_WORLD, &request);
                       MPI_Wait(&request, &status);
                       //printf("rank: %d, LOOP3 step 1-I\n",my_rank);
                       //fflush(stdout);
                       goto LOOP2;
                   }
                } else {
                                           // if no incoming job requests
                                           // start requesting
                   //printf("rank: %d sending message to %d, LOOP3 step 1-1\n",my_rank,dest);
                   //fflush(stdout);
                   MPI_Isend(&my_rank,1,MPI_INT,dest,4,MPI_COMM_WORLD,&request);
                   MPI_Wait(&request,&status);
                                           // request a buddy and busy-wait
                                           // on its reply
                  while (1) {
                   //printf("rank: %d, LOOP3 step 1-B, waiting on %d\n",my_rank, dest);
                   //fflush(stdout);
                                           // if received a job request from
                                           // any other buddy, send neg reply
                   flag = 0;
                   MPI_Iprobe(MPI_ANY_SOURCE,4,MPI_COMM_WORLD,&flag,&status);
                   //printf("rank: %d, LOOP3 step 1-1-A, waiting on %d\n",my_rank, dest);
                   //fflush(stdout);
                   if (flag) {
                   //printf("rank: %d, LOOP3 step 1-1-B, recd job reqs from %d\n",my_rank, nBC);
                   //fflush(stdout);
                    MPI_Irecv(&nBC, 1, MPI_INT, MPI_ANY_SOURCE,                \
                              4, MPI_COMM_WORLD, &request);
                    MPI_Wait(&request, &status);
                    MPI_Isend(&my_rank,1,MPI_INT,nBC,6,MPI_COMM_WORLD,&request);
                    MPI_Wait(&request,&status);
                   //printf("rank: %d, LOOP3 step 1-1-C, sending job na to %d\n",my_rank, nBC);
                   //fflush(stdout);
                   }
                                           // if no job requests pending, check
                                           // for positive reply
                   flag = 0;
                   MPI_Iprobe(dest,5,MPI_COMM_WORLD,&flag,&status);
                   //printf("rank: %d, LOOP3 step 1-C\n",my_rank);
                   //fflush(stdout);
                   if (flag==0) {
                                           // if no positive reply, check for
                                           // negative reply
                    //printf("rank: %d, LOOP3 step 1-D probing -r frm %d\n",my_rank, dest);
                    //fflush(stdout);
                    MPI_Iprobe(dest,6,MPI_COMM_WORLD,&flag,&status);
                    //printf("rank: %d, LOOP3 step 1-E probing -r frm %d\n",my_rank,dest);
                    //fflush(stdout);
                    if (flag==0) continue;

                                           // if negative reply, receive it, 
                                           // incr mnr, & move to next buddy
                    //printf("rank: %d, LOOP3 step 1-F\n",my_rank);
                    //fflush(stdout);
                    MPI_Irecv(&nBC, 1, MPI_INT, dest,                         \
                              6, MPI_COMM_WORLD, &request);
                    MPI_Wait(&request, &status);
                    //printf("rank: %d, LOOP3 step 1-G\n",my_rank);
                    //fflush(stdout);
                    ++mnr;
                                           // while moving to next buddy loop
                                           // around if required
                    ++dest;
                    if (dest==numProcs-1) {
                       dest=0;
                    }
                                           // if no job in system, wait to 
                                           // terminate
                    if (mnr==numProcs-1) {
                       //printf("rank: %d, LOOP3 step 1-H\n",my_rank);
                       //fflush(stdout);
                       MPI_Isend(&my_rank, 1, MPI_INT, numProcs-1,            \
                                 2, MPI_COMM_WORLD, &request);
                       MPI_Wait(&request, &status);
                       //printf("rank: %d, LOOP3 step 1-I\n",my_rank);
                       //fflush(stdout);
                       goto LOOP2;
                    }
                    break;
                   } else {
                                           // if positive reply to job request
                                           // receive, unpack and process
                    //printf("rank: %d, LOOP3 step 1-J\n",my_rank);
                    //fflush(stdout);
                    MPI_Irecv(job, sopam, MPI_INT, dest,                     \
                              5, MPI_COMM_WORLD, &request);
                    MPI_Wait(&request, &status);
                    //printf("rank: %d, LOOP3 step 1-K\n",my_rank);
                    //fflush(stdout);
                    my_level = job[0];
                    //printf("rank: %d, LOOP3 step 1-K-1\n",my_rank);
                    //fflush(stdout);
                    lsize = job[1];
                    //printf("rank: %d, LOOP3 step 1-K-2\n",my_rank);
                    //fflush(stdout);
                    rsize = job[2];
                    //printf("rank: %d, LOOP3 step 1-K-3\n",my_rank);
                    //fflush(stdout);
                    memcpy(workSpace[my_level][0][0],&job[1],
                                                            lsize*sizeof(int));
                    //printf("rank: %d, LOOP3 step 1-K-4\n",my_rank);
                    //fflush(stdout);
                    memcpy(workSpace[my_level][0][1],&job[1+lsize],
                                                            rsize*sizeof(int));
                    //printf("rank: %d, LOOP3 step 1-L\n",my_rank);
                    //fflush(stdout);

                    nctr = 0;
                    taxa = taxaQueue[job[6]];
                    //printf("rank: %d, LOOP3 step 1-M, my level: %d\n",my_rank,my_level);
                    //fflush(stdout);

                    for (indx=lhsOffSet+skipSize; indx<=lsize-skipSize;       \
                                                            indx+=skipSize) {
                          genTaskLeft(workSpace[my_level][0],                 \
                                        workSpace[my_level+1][nctr],taxa,indx);
                          nctr = nctr + 1;
                    }
                    //printf("rank: %d, LOOP3 step 1-N\n",my_rank);
                    //fflush(stdout);
                    for (indx=1+skipSize;indx<=rsize-skipSize;indx+=skipSize) {
                          genTaskRight(workSpace[my_level][0],                \
                                        workSpace[my_level+1][nctr],taxa,indx);
                          nctr = nctr + 1;
                    }
                    wosac[my_level+1][0] = 0;
                    //printf("rank: %d, LOOP3 step 1-O\n",my_rank);
                    //fflush(stdout);
                    my_level = my_level+1;
                    goto LOOP1;
                   }
                  }
                }
             }
           }
           //printf("rank: %d, LOOP3 step 2\n",my_rank);
           //fflush(stdout);
  }

TERMINATE:
  numProcs = numProcs;
}
