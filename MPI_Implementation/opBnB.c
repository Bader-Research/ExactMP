#include "header.h"

void paBnB(void) {
  int level=0, rsize, lsize, cctr, nctr, taxa, cost, indx, mnr=0;
  
  int flag, dest, nbr, oBC, nBC;
                                           // If I am the last processor, I
                                           // wouldnt do crap work :)
  if (my_rank==numProcs-1) {
    while (1) {
      flag = 0;
      MPI_Iprobe(MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&flag,&status);
      if (flag==1) {
         MPI_Irecv(&nBC,1,MPI_INT,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&request);
         MPI_Wait(&request, &status);
         ++mnr;
      }
      if (mnr==numProcs-1) {
         for (cctr=0; cctr < numProcs-1; cctr++) {
            MPI_Isend(&my_rank,1,MPI_INT,cctr,1,MPI_COMM_WORLD,&request);
            MPI_Wait(&request, &status);
         }
         goto TERMINATE;
      }
    }
  } else {
                                           // Else, I do the main job of   
                                           // computing the global minima
    initBestCost = bestCost;
    numBestSolutions=0;
    //bestCost = 600;
                                           // Only if I have a job at hand,
                                           // will I get in...
    if (my_level >= 0) {
LOOP1: 
      level = my_level;
      mnr = 0;                             // my negative replies = 0 as of now

                                           // as long as I have jobs at hand, I
                                           // slog on...
      while (wosac[my_level][0] < wosac[my_level][1]) {
                                           // record current best cost to 
                                           // propagate if it improves
        oBC = bestCost;
                                           // if this is the last level,  
                                           // check for best cost change
        if (level==matrix.num_taxa-7) {
          for (indx=0; indx < wosac[level][1]; indx++) {
            if (workSpace[level][indx][0][8] < bestCost) {
              bestCost = workSpace[level][indx][0][8];
	      printf("\t\t\t\t\t\tRank[%d]: [BEST SCORE: %d]\n",my_rank,      \
                                                             bestCost+pNiCost);
              fflush(stdout);
                                                   // Copy LHS Tree
              //lsize = workSpace[level][indx][0][0];
              //memcpy(solQ[0][0],workSpace[level][indx][0],lsize*sizeof(int));
                                                   // Copy RHS tree
              //rsize = workSpace[level][indx][0][1];
              //memcpy(solQ[0][1],workSpace[level][indx][1],rsize*sizeof(int));
      	      numBestSolutions = 1;
	    } 
          }
          wosac[level][0]=wosac[level][1];
          if (level > my_level) {
	    --level;
          }
        } else {
                                           // if this is not last level,  
                                           // expand job into next level & loop
          cctr = wosac[level][0];
          lsize = workSpace[level][cctr][0][0];
          rsize = workSpace[level][cctr][0][1];
          nctr = 0;
          taxa = taxaQueue[workSpace[level][cctr][0][5]];
          //cost = workSpace[level][cctr][0][8];
          cost = workSpace[level][cctr][0][8] + plbValues[level];
          if (cost < bestCost) {
            for (indx=lhsOffSet+skipSize;indx<=lsize-skipSize;indx+=skipSize){
              genTaskLeft(workSpace[level][cctr],                             \
                                           workSpace[level+1][nctr],taxa,indx);
              nctr = nctr + 1;
            }
            for (indx=1+skipSize;indx<=rsize-skipSize;indx+=skipSize) {
              genTaskRight(workSpace[level][cctr],                            \
                                           workSpace[level+1][nctr],taxa,indx);
              nctr = nctr + 1;
            }
            wosac[level+1][0] = 0;
            ++wosac[level][0];
            ++level;
          } else {
             ++wosac[level][0];
          }
        }
                                           // check if this job was last job
        for (cctr=level; cctr > my_level; cctr--) {
          if (wosac[level][0]==wosac[level][1]) {
	    --level;
          } else {
            break;
          }
        }

                    // PARALLEL STUFF: check messages for me
                                           // if better best cost, propagate it
                                           // to all buddies, but not to self
        if (bestCost < oBC) {
          for (dest=0; dest < numProcs-1; dest++) {
            if (dest!=my_rank) {
              MPI_Isend(&bestCost,1,MPI_INT,dest,3,MPI_COMM_WORLD,&request);
              MPI_Wait(&request, &status);
            }
          }
        }
                                           // if any incoming messages are here
                                           // finish them, before proceeding
        while (1) {
                                           // tag=3 : better-bestcost message
          flag = 0;
          MPI_Iprobe(MPI_ANY_SOURCE,3,MPI_COMM_WORLD,&flag,&status);
          if (flag==1) {
            MPI_Irecv(&nBC,1,MPI_INT,MPI_ANY_SOURCE,3,MPI_COMM_WORLD,&request);
            MPI_Wait(&request, &status);
            if (nBC < bestCost) {
              bestCost = nBC;
            }
          }
                                           // tag=4 : job request message
          flag = 0;
          MPI_Iprobe(MPI_ANY_SOURCE,4,MPI_COMM_WORLD,&flag,&status);
          if (flag==1) {
           MPI_Irecv(&dest,1,MPI_INT,MPI_ANY_SOURCE,4,MPI_COMM_WORLD,&request);
           MPI_Wait(&request, &status);
                                           // if jobs at hand, pack and send
           if (wosac[my_level][0]<wosac[my_level][1] &&                       \
                                                  my_level<matrix.num_taxa-7) {
             indx = wosac[my_level][0];
             lsize = workSpace[my_level][indx][0][0];
             rsize = workSpace[my_level][indx][0][1];
             job[0] = my_level;
             memcpy(&job[1],workSpace[my_level][indx][0],lsize*sizeof(int));
             memcpy(&job[1+lsize],workSpace[my_level][indx][1],               \
                                                            rsize*sizeof(int));
                                           // tag=5 : job available and sending
             MPI_Isend(job,sopam,MPI_INT,dest,5,MPI_COMM_WORLD,&request);
             MPI_Wait(&request, &status);
             ++wosac[my_level][0];
           } else {
                                           // if no jobs at hand send neg reply
                                           // tag=6
             MPI_Isend(&my_rank,1,MPI_INT,dest,6,MPI_COMM_WORLD,&request);
             MPI_Wait(&request, &status);
           }
          } else if (flag==0) {
                                           // if no messages, move on...
             break;
          }
        } // WHILE (1)
      } // UNTIL HERE IT WAS BNB PROCESSING
    } // if (my_level >= 0)

                    // PARALLEL STUFF: loop here after retirement
    while (1) {
LOOP2:
                                           // probe & act on terminate message
      flag = 0;
      MPI_Iprobe(numProcs-1,1,MPI_COMM_WORLD,&flag,&status);
      if (flag==1) {
        MPI_Irecv(&dest,1,MPI_INT,numProcs-1,1,MPI_COMM_WORLD,&request);
        MPI_Wait(&request, &status);
        goto TERMINATE;
      }
                                           // probe & reply job request message
      flag = 0;
      MPI_Iprobe(MPI_ANY_SOURCE,4,MPI_COMM_WORLD,&flag,&status);
      if (flag==1) {
        MPI_Irecv(&dest,1,MPI_INT,MPI_ANY_SOURCE,4,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,&status);
        MPI_Isend(&my_rank,1,MPI_INT,dest,6,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,&status);
      }
                                           // probe & act bestcost message
      flag = 0;
      MPI_Iprobe(MPI_ANY_SOURCE,3,MPI_COMM_WORLD,&flag,&status);
      if (flag==1) {
        MPI_Irecv(&nBC,1,MPI_INT,MPI_ANY_SOURCE,3,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,&status);
        if (nBC < bestCost) {
          bestCost = nBC;
        }
      }

                                           // if all buddies done, wait on
                                           // terminate message from slughead
      if (mnr==numProcs-1) {
        goto LOOP2; 
      }
                                           // if no messages, and more jobs
                                           // available, move on...
      if (flag==0) {
        break;
      }
   }

                    // PARALLEL STUFF: request a job, if none wait to terminate
LOOP3:
                                           // goto neighbor first, loop if reqd
   if ((dest=my_rank+1)==numProcs-1) dest=0;
                                           // request all buddies for job
   for (mnr=1, nbr=0; nbr < numProcs-1; nbr++) {
    if (nbr!=my_rank) {
     flag = 0;
                                           // first check for any incoming 
                                           // job requests & send neg replies
     MPI_Iprobe(dest,4,MPI_COMM_WORLD,&flag,&status);
     if (flag==1) {
      MPI_Irecv(&nBC,1,MPI_INT,dest,4,MPI_COMM_WORLD,&request);
      MPI_Wait(&request, &status);
      MPI_Isend(&my_rank,1,MPI_INT,dest,6,MPI_COMM_WORLD,&request);
      MPI_Wait(&request,&status);
      ++mnr;
                                           // if no job in system, wait to 
                                           // terminate
      if (mnr==numProcs-1) {
       MPI_Isend(&my_rank,1,MPI_INT,numProcs-1,2,MPI_COMM_WORLD,&request);
       MPI_Wait(&request, &status);
       goto LOOP2;
      }
     } else {
                                           // if no incoming job requests
                                           // start requesting
       MPI_Isend(&my_rank,1,MPI_INT,dest,4,MPI_COMM_WORLD,&request);
       MPI_Wait(&request,&status);
                                           // request a buddy and busy-wait
                                           // on its reply
       while (1) {
                                           // if received a job request from
                                           // any other buddy, send neg reply
         flag = 0;
         MPI_Iprobe(MPI_ANY_SOURCE,4,MPI_COMM_WORLD,&flag,&status);
         if (flag) {
           MPI_Irecv(&nBC,1,MPI_INT,MPI_ANY_SOURCE,4,MPI_COMM_WORLD,&request);
           MPI_Wait(&request, &status);
           MPI_Isend(&my_rank,1,MPI_INT,nBC,6,MPI_COMM_WORLD,&request);
           MPI_Wait(&request,&status);
         }
                                           // if no job requests pending, check
                                           // for positive reply
         flag = 0;
         MPI_Iprobe(dest,5,MPI_COMM_WORLD,&flag,&status);
         if (flag==0) {
                                           // if no positive reply, check for
                                           // negative reply
           MPI_Iprobe(dest,6,MPI_COMM_WORLD,&flag,&status);
           if (flag==0) continue;

                                           // if negative reply, receive it, 
                                           // incr mnr, & move to next buddy
           MPI_Irecv(&nBC,1,MPI_INT,dest,6,MPI_COMM_WORLD,&request);
           MPI_Wait(&request, &status);
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
            MPI_Isend(&my_rank,1,MPI_INT,numProcs-1,2,MPI_COMM_WORLD,&request);
            MPI_Wait(&request, &status);
            goto LOOP2;
           }
           break;
         } else {
                                           // if positive reply to job request
                                           // receive, unpack and process
           MPI_Irecv(job,sopam,MPI_INT,dest,5,MPI_COMM_WORLD,&request);
           MPI_Wait(&request, &status);
           my_level = job[0];
           lsize = job[1];
           rsize = job[2];
           memcpy(workSpace[my_level][0][0],&job[1],lsize*sizeof(int));
           memcpy(workSpace[my_level][0][1],&job[1+lsize],rsize*sizeof(int));

           nctr = 0;
           taxa = taxaQueue[job[6]];

           for (indx=lhsOffSet+skipSize;indx<=lsize-skipSize;indx+=skipSize) {
             genTaskLeft(workSpace[my_level][0],                              \
                                        workSpace[my_level+1][nctr],taxa,indx);
             nctr = nctr + 1;
           }
           for (indx=1+skipSize;indx<=rsize-skipSize;indx+=skipSize) {
             genTaskRight(workSpace[my_level][0],                             \
                                        workSpace[my_level+1][nctr],taxa,indx);
             nctr = nctr + 1;
           }
           wosac[my_level+1][0] = 0;
           my_level = my_level+1;
           goto LOOP1;
         }
       }
     }
    }
   }
  }

TERMINATE:
  numProcs = numProcs;
}
