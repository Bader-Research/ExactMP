#include "header.h"

struct maxMiniTaxaPosWLength getEckDayMinTree(struct subTask *sourceTree, int taxaNum) {
	struct subTask *childTree;
	int currentLength = 0;
	int prevLength = 0;
 
	struct maxMiniTaxaPosWLength mMTPWL; 
	 
	int pos;
	int currMinPos; 
	int currMinTaxa; 
	int arraySize = sourceTree->sizeOfSubTaskArray; 
	 
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
	printf("\n\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	printf("\t\t\t\t CHECKING MINIMUM TREE LENGTH FOR ALL POSITIONS FOR TAXA: %d\n",taxaNum);
	printf("\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
#endif
                                                             /* generate first case, adding taxa in position 1 */
	childTree = genNextEckDayTask(sourceTree, childTree, taxaNum, 1); 
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
	printTree(childTree, 0);
#endif
                                                             /* update initial tree length to the length in position 1 */
	prevLength = childTree->costOfSubTask;
	freeTree(childTree);

	currMinPos = 1;
	currMinTaxa = taxaNum; 

                                                             /* for all remaining positions */
	for (pos=2; pos < arraySize; pos++) { 
                                                             /* do not try position before internal node -1 */
                                                             /* this is done to ensure that the integrity of tree */
			childTree = genNextEckDayTask(sourceTree,childTree, taxaNum, pos); 
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
	                printTree(childTree, pos);
#endif
			currentLength = childTree->costOfSubTask; 
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
			printf("\n\t\t->->->->->->RETAIN MINIMUM LENGTH FROM CURRENT & PREVIOUS LENGTHS<-<-<-<-<-<-\n");
                        printf("\n\t\t\tCurr {MinLen:%3d @ Pos:%3d w/ Taxa:%3d} vs. Prev {MinLen:%3d @ Pos:%3d w/ Taxa:%3d}\n\n",currentLength, pos, taxaNum, prevLength, currMinPos, currMinTaxa);
#endif
                                                             /* if this length is less than previous length, keep it */
			if (currentLength < prevLength) {
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
				printf("\t\tCurr Lgt (%d) < Prev Lgt (%d)=> ",currentLength,prevLength);
#endif
				prevLength = currentLength; 
				currMinPos = pos; 
				currMinTaxa = taxaNum; 
#ifdef ECK_AND_DAY_MIN_TREE_VERBOSE
				printf("Update Curr Lgt to: %d & keeping Taxa: %d, Pos: %d\n\n",currentLength,currMinTaxa,currMinPos);
#endif
			} 

			freeTree(childTree);
	}
		 
                                                             /* return position & tree with minimum length */
	mMTPWL.taxa = currMinTaxa;
	mMTPWL.pos = currMinPos;
	mMTPWL.length = prevLength;
		 
	return (mMTPWL); 
}

void computeEckDayCost (struct subTask *sourceSubTask) {

  int rChildIndx, lChildIndx, arrIndx;
  struct subTaskArrayNode *leftChild, *rightChild, *parent;

                                                         /* Loop until root of the subtree is encountered */
  for (arrIndx=sourceSubTask->sizeOfSubTaskArray-1; arrIndx >= 1; arrIndx--) {
                                                         /* For all internal nodes which are also parent nodes, */
	if (sourceSubTask->node[arrIndx].id < 0) {
                                                         /* Gather children information */
		lChildIndx = sourceSubTask->node[arrIndx].lChildIndx;
		leftChild = &sourceSubTask->node[lChildIndx];

		rChildIndx = sourceSubTask->node[arrIndx].rChildIndx;
		rightChild = &sourceSubTask->node[rChildIndx];

		parent = &sourceSubTask->node[arrIndx];
                                                         /* Send left child, parent, right child for Fitch's operation */
		performFitchOp(leftChild , parent, rightChild);
	}
  }
                                                         /* Both subtrees (left & right) are computed and cost kept */
                                                         /* in the roots of subtrees, compute cost of root node now */
  sourceSubTask->root.rChildIndx = 1;
  sourceSubTask->root.lChildIndx = 0;
  performFitchOp(&sourceSubTask->node[0] , &sourceSubTask->root, &sourceSubTask->node[1]);
  sourceSubTask->costOfSubTask = sourceSubTask->root.costUntilThisNode;
}

void eckDay(void) {
	int taxaOver = 3, ctr, stackPtr, lChildIndx, rChildIndx, *locationStack, keepLoc;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;

	sourceTree = genIstSubTask();
                                                             /* apply maxmini, & enqueue first addition order taxa */
	taxaOver += 1; 
	mMTPWL = getEckDayMinTree(sourceTree, taxaOver);
                                                             /* for all taxa remaining to be enqueued in addition order */
	while (taxaOver < matrix.num_taxa) {
#ifdef ECK_AND_DAY_VERBOSE
		printf("\t\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~\n");
		printf("\t\t\t\t\t\tADDING NEXT TAXA ORDER\n");
		printf("\t\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~");
#endif
                                                             /* generate next task with previous maxmini taxa, position */
		childTree = genNextEckDayTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos); 
                                                             /* free parent task, it is not needed anymode */
		freeTree(sourceTree);
		sourceTree = childTree;
#ifdef ECK_AND_DAY_VERBOSE
		printTree(sourceTree,mMTPWL.pos);
#endif
                                                             /* get next position of minimum cost */
		taxaOver += 1; 
		mMTPWL = getEckDayMinTree(sourceTree, taxaOver);
	}
                                                             /* if the final tree cost is less than best cost computed */
  	childTree = genNextEckDayTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos); 

	if (childTree->costOfSubTask < bestCost) {
	                                                     /* The following lines of code just make sure that */
                                                             /* the root of RHS sub tree has an internal ID of -1 */
	                                                     /* this is done because, other functions assume -1 to be */
	                                                     /* root of RHS tree. The line goes from here until */
		stackPtr = -1;
		locationStack = (int *)calloc(childTree->sizeOfSubTaskArray, sizeof(int));
		for (ctr=childTree->sizeOfSubTaskArray-1; ctr >= 0; ctr--) {
			if (childTree->node[ctr].id > 0) {
				stackPtr = stackPtr + 1;
				locationStack[stackPtr] = ctr;
			}
		                                        // if id < 0, copy state vectors and pop twice, compute cost and
							// push the negative node onto stack
			else if (childTree->node[ctr].id < 0) {
				if (childTree->node[ctr].id == -1) {
					keepLoc = ctr;
				}
				lChildIndx = locationStack[stackPtr];
				stackPtr = stackPtr - 1;
	
				rChildIndx = locationStack[stackPtr];
				stackPtr = stackPtr - 1;
	
				stackPtr = stackPtr + 1;
				locationStack[stackPtr] = ctr;
			}
		  }
		  lChildIndx = locationStack[stackPtr];
		  stackPtr = stackPtr - 1;

		  rChildIndx = locationStack[stackPtr];
		  stackPtr = stackPtr - 1;
	
		  if (childTree->node[rChildIndx].id != -1) {
			ctr = childTree->node[rChildIndx].id;
			childTree->node[rChildIndx].id = childTree->node[keepLoc].id;
			childTree->node[keepLoc].id = ctr;
		  } 
	                                                     /* setting RHS subtree root id to be -1 ends here */
		  free(locationStack);
	                                                     /* by ED algorithm, set this as best cost */
		  bestCost = childTree->costOfSubTask;
		  flushBestCostStack();
		  freeTree(sourceTree);
		  pushBestCostNode(childTree);
	}
	else {
		  freeTree(childTree);
		  freeTree(sourceTree);
	}
}

struct subTask *genNextEckDayTask(struct subTask *currTask, struct subTask *nextTask, int taxa, int pos) {

   int i;

   nextTask = (struct subTask *)calloc(1, sizeof(struct subTask));

                                                        /* get highest internal node of the child task */
   nextTask->rhtSubTreeRootIndx = 1;
   nextTask->highestInternalNode = currTask->highestInternalNode-1;
                                                        /* get next taxa to be added */
   nextTask->lastTaxaAdded = taxa;
                                                        /* get # elements in the preordered array of parent task */
   nextTask->sizeOfSubTaskArray = currTask->sizeOfSubTaskArray+2;

   nextTask->node = (struct subTaskArrayNode *)calloc(nextTask->sizeOfSubTaskArray,sizeof(struct subTaskArrayNode));

                                                        /* copy all elements from 0 to pos from parent to child task */
   for (i = 0; i < pos; i++) {
     nextTask->node[i].id = currTask->node[i].id; 
     nextTask->node[i].costUntilThisNode = currTask->node[i].costUntilThisNode; 
     nextTask->node[i].costAtThisNode = currTask->node[i].costAtThisNode; 
     nextTask->node[i].stateVector = (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
     memcpy(nextTask->node[i].stateVector, currTask->node[i].stateVector, matrix.num_pars_inf_sites*sizeof(int));
   }

                                                        /* 2 elems from "pos" get new values (internal node, taxa) pair */
   nextTask->node[pos].id = nextTask->highestInternalNode;
   nextTask->node[pos].costUntilThisNode = 0;
   nextTask->node[pos].costAtThisNode = 0;
   nextTask->node[pos].stateVector = (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));

   nextTask->node[pos+1].id = nextTask->lastTaxaAdded;
   nextTask->node[pos+1].costUntilThisNode = 0;
   nextTask->node[pos+1].costAtThisNode = 0;
   nextTask->node[pos+1].stateVector = (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
   memcpy(nextTask->node[pos+1].stateVector, matrix.reord_sites_enc[nextTask->lastTaxaAdded-1], matrix.num_pars_inf_sites*sizeof(int));
                                                        /* copy all elements from "pos" in parent to "pos+2" in child */
   for (i=pos; i < currTask->sizeOfSubTaskArray; i++) {
     nextTask->node[i+2].id = currTask->node[i].id; 
     nextTask->node[i+2].costUntilThisNode = currTask->node[i].costUntilThisNode; 
     nextTask->node[i+2].costAtThisNode = currTask->node[i].costAtThisNode; 
     nextTask->node[i+2].stateVector = (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
     memcpy(nextTask->node[i+2].stateVector, currTask->node[i].stateVector, matrix.num_pars_inf_sites*sizeof(int));
   }

   nextTask->root.stateVector = (int *)calloc(matrix.num_pars_inf_sites, sizeof(int));
                                                        /* set parent child relationship for the child tree */
   setParentChildIndices(nextTask);
                                                        /* compute cost of the new child */
   computeEckDayCost(nextTask);
   return (nextTask);

}


