#include "header.h"

                                                              /* This function employs 3 randomized methods to search */
                                                              /* for better incumbent. This is deployed after a best */
                                                              /* cost is found by Eck-Day, Modified Greedy, NJ, and  */
                                                              /* Maxmini algorithms */
void randOptimize(void) {
	int ctr=0, **arr, numTbrTrees, randNum, oldBc;

	                                                      /* generate arbitrary addition order sequences, number */
	                                                      /* of such sequences can be specified in header file */
	arr = (int **)calloc(NUM_TBR_TREES, sizeof(int *));
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		arr[numTbrTrees] = (int *)calloc(matrix.num_taxa-3, sizeof(int));
		for (ctr=0; ctr < matrix.num_taxa-3; ctr++) {
			arr[numTbrTrees][ctr] = -1;
		}
	}

	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		ctr=0;
		//srandom(numTbrTrees+1);
		srandom(random()+1);
		while (ctr < matrix.num_taxa-3) {
			randNum	= random()%(matrix.num_taxa-3);
			if (arr[numTbrTrees][randNum] == -1) {
				arr[numTbrTrees][randNum] = ctr+4;
				++ctr;
			}
		}
	}

	                                                      /* apply eckDay algorithm to each of those addition */
	                                                      /* orders for better solution */
    if (randOptLev==1) {
	oldBc = bestCost+pNiCost;
	printf("\tRandomization level I:");
	fflush(stdout);
	//fprintf(fp,"\tRandomization level I:");
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		nQ4RandOpt(arr, numTbrTrees);
		TBRSmartGreedy();
		printf("\n\t\t[Iteration-> %3d ]",numTbrTrees);
		fflush(stdout);
		if (oldBc > bestCost+pNiCost) {
			printf("\t\t[BEST SCORE %d]",bestCost+pNiCost);
			fflush(stdout);
			//fprintf(fp,"\t\t[BEST SCORE %d]",bestCost+pNiCost);
			oldBc = bestCost+pNiCost;
		}
	}
	printf("\n");
	fflush(stdout);
    }

	                                                      /* apply eckDay, and maxmini algorithms to each of those */
	                                                      /* addition orders for better solution */
    if (randOptLev==2) {
	oldBc = bestCost+pNiCost;
	printf("\tRandomization level I:");
	fflush(stdout);
	//fprintf(fp,"\tRandomization level I:");
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		nQ4RandOpt(arr, numTbrTrees);
		TBRSmartGreedy();
		printf("\n\t\t[Iteration-> %3d ]",numTbrTrees);
		fflush(stdout);
		if (oldBc > bestCost+pNiCost) {
			printf("\t\t[BEST SCORE %d]",bestCost+pNiCost);
			fflush(stdout);
			//fprintf(fp,"\t\t[BEST SCORE %d]",bestCost+pNiCost);
			oldBc = bestCost+pNiCost;
		}
	}
	printf("\n");
	fflush(stdout);

	oldBc = bestCost+pNiCost;
	printf("\tRandomization level II:");
	fflush(stdout);
	//fprintf(fp,"\tRandomization level II:");
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		TBReckDay(arr, numTbrTrees);
		printf("\n\t\t[Iteration-> %3d ]",numTbrTrees);
		fflush(stdout);
		if (oldBc > bestCost+pNiCost) {
			printf("\t\t[BEST SCORE %d]",bestCost+pNiCost);
			fflush(stdout);
			//fprintf(fp,"\t\t[BEST SCORE %d]",bestCost+pNiCost);
			oldBc = bestCost+pNiCost;
		}
	}
	printf("\n");
	fflush(stdout);
    }

	                                                      /* apply eckDay, maxmini, and modified greedy algorithms */
	                                                      /* each of those addition orders for better solution */
    if (randOptLev==3) {
	oldBc = bestCost+pNiCost;
	printf("\tRandomization level I:");
	fflush(stdout);
	//fprintf(fp,"\tRandomization level I:");
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		nQ4RandOpt(arr, numTbrTrees);
		TBRSmartGreedy();
		printf("\n\t\t[Iteration-> %3d ]",numTbrTrees);
		fflush(stdout);
		if (oldBc > bestCost+pNiCost) {
			printf("\t\t[BEST SCORE %d]",bestCost+pNiCost);
			fflush(stdout);
			//fprintf(fp,"\t\t[BEST SCORE %d]",bestCost+pNiCost);
			oldBc = bestCost+pNiCost;
		}
	}
	printf("\n");
	fflush(stdout);

	oldBc = bestCost+pNiCost;
	printf("\tRandomization level II:");
	fflush(stdout);
	//fprintf(fp,"\tRandomization level II:");
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		TBReckDay(arr, numTbrTrees);
		printf("\n\t\t[Iteration-> %3d ]",numTbrTrees);
		fflush(stdout);
		if (oldBc > bestCost+pNiCost) {
			printf("\t\t[BEST SCORE %d]",bestCost+pNiCost);
			fflush(stdout);
			//fprintf(fp,"\t\t[BEST SCORE %d]",bestCost+pNiCost);
			oldBc = bestCost+pNiCost;
		}
	}
	printf("\n");
	fflush(stdout);

	oldBc = bestCost+pNiCost;
	printf("\tRandomization level III:");
	fflush(stdout);
	//fprintf(fp,"\tRandomization level III:");
	for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
		nQ4RandOpt(arr, numTbrTrees);
		TBRmaxMin();
		printf("\n\t\t[Iteration-> %3d ]",numTbrTrees);
		fflush(stdout);
		if (oldBc > bestCost+pNiCost) {
			printf("\t\t[BEST SCORE %d]",bestCost+pNiCost);
			fflush(stdout);
			//fprintf(fp,"\t\t[BEST SCORE %d]",bestCost+pNiCost);
			oldBc = bestCost+pNiCost;
		}
	}
	printf("\n");
	fflush(stdout);
    }

    for (numTbrTrees=0; numTbrTrees < NUM_TBR_TREES; numTbrTrees++) {
	free(arr[numTbrTrees]);
    }
    free(arr);

}

void TBReckDay(int **arr, int row) {
	int taxaOver = 3, ctr;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;

	sourceTree = genIstSubTask();
                                                             /* apply maxmini, & enqueue first addition order taxa */
	ctr = 0;
	taxaOver += 1; 
	mMTPWL = getEckDayMinTree(sourceTree, arr[row][ctr]);
	ctr = ctr + 1;
                                                             /* for all taxa remaining to be enqueued in addition order */
	while (taxaOver < matrix.num_taxa) {
                                                             /* generate next task with previous maxmini taxa, position */
		childTree = genNextEckDayTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos); 
                                                             /* free parent task, it is not needed anymode */
		freeTree(sourceTree);
		sourceTree = childTree;
                                                             /* get next position of minimum cost */
		taxaOver += 1; 
		mMTPWL = getEckDayMinTree(sourceTree, arr[row][ctr]);
		ctr = ctr + 1;
	}
                                                             /* if the final tree cost is less than best cost computed */
  	childTree = genNextEckDayTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos); 
	if (childTree->costOfSubTask < bestCost) {
		  bestCost = childTree->costOfSubTask;
	}

	correctEckDayIndices(childTree);
	tbr(childTree);

	freeTree(childTree);
	freeTree(sourceTree);
}

void correctEckDayIndices (struct subTask *sourceSubTask) {
  int arrIndx, indx, indxOfMinusOne;


                                                        /* for all nodes in the array starting from RHS */
  for (arrIndx = sourceSubTask->sizeOfSubTaskArray-1; arrIndx >= 0; arrIndx = arrIndx-1) {
                                                        /* if node is a leaf node, push it to stack */
    if (sourceSubTask->node[arrIndx].id > 0) {
      sourceSubTask->node[arrIndx].rChildIndx = LEAF;
      sourceSubTask->node[arrIndx].lChildIndx = LEAF;
      push(arrIndx);
    }
                                                        /* if node is internal node, pop 2 nodes from stack, & push this */
                                                        /* node on stack. 2 poped nodes are left & right children resp. */
    else if (sourceSubTask->node[arrIndx].id < 0) {
      if (sourceSubTask->node[arrIndx].id == -1) {
        indxOfMinusOne = arrIndx;
      }

      indx = pop();
      sourceSubTask->node[arrIndx].lChildIndx = indx;
      sourceSubTask->node[indx].parentIndx = arrIndx;

      indx = pop();
      sourceSubTask->node[arrIndx].rChildIndx = indx;
      sourceSubTask->node[indx].parentIndx = arrIndx;

      push(arrIndx);
    }
  }
                                                        /* Children of main root are roots of rhs & lhs subtrees */
  sourceSubTask->root.lChildIndx = pop();
  sourceSubTask->root.rChildIndx = pop();

  indx = sourceSubTask->node[sourceSubTask->root.rChildIndx].id;

  sourceSubTask->node[sourceSubTask->root.rChildIndx].id = -1;
  sourceSubTask->node[indxOfMinusOne].id = indx;
  sourceSubTask->rhtSubTreeRootIndx = sourceSubTask->root.rChildIndx;
                                                        /* Parent of roots of rhs & lhs subtrees is main root */
  sourceSubTask->node[0].parentIndx = ROOT;
  sourceSubTask->node[sourceSubTask->rhtSubTreeRootIndx].parentIndx = ROOT;

}

void nQ4RandOpt(int **arr, int row) {
  struct taxaQ* temp_node;
  int ctr;
  for (ctr=0; ctr < matrix.num_taxa-3; ctr++) {
    if (taxaQHead==NULL) {
        taxaQHead = (struct taxaQ*)malloc(sizeof(struct taxaQ));
        taxaQHead->next = NULL;
      	taxaQHead->taxa = arr[row][ctr];
        taxaQTail = taxaQHead;
    }
    else {
        temp_node = (struct taxaQ*)malloc(sizeof(struct taxaQ));
        taxaQTail->next = temp_node;
      	temp_node->taxa = arr[row][ctr];
        temp_node->next = NULL;
        taxaQTail = temp_node;
    }
  }
}

void TBRmaxMin(void) {
	int taxaOver = 3;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;

	sourceTree = genIstSubTask();
                                                             /* apply maxmini, & enqueue first addition order taxa */
	mMTPWL = maxMini(sourceTree);
	taxaOver += 1;
                                                             /* for all taxa remaining to be enqueued in addition order */
	while (taxaOver < matrix.num_taxa) {
                                                             /* generate next task with previous maxmini taxa, position */
		childTree = genNextTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos);
                                                             /* free parent task, it is not needed anymode */
		freeTree(sourceTree);
		sourceTree = childTree;
                                                             /* apply maxmini, & enqueue next addition order taxa */
		mMTPWL = maxMini(sourceTree);
		taxaOver += 1;
	}

	childTree = genNextTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos);
	freeTree(sourceTree);
	tbr(childTree);
                                                             /* if the final tree cost is less than best cost computed */
	                                                     /* by SG algorithm, set this as best cost */
	if (childTree->costOfSubTask < bestCost) {
	  bestCost = childTree->costOfSubTask;
	}
	freeTree(childTree);
}

void TBRSmartGreedy(void) {
	int taxaOver = 3;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;

	sourceTree = genIstSubTask();
	mMTPWL = smartGreedy(sourceTree); 
	taxaOver += 1; 

	while (taxaOver < matrix.num_taxa) { 
		childTree = genNextTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos); 
		freeTree(sourceTree);
		sourceTree = childTree;
		mMTPWL = smartGreedy(sourceTree);
		taxaOver += 1; 
	}

	childTree = genNextTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos);
	freeTree(sourceTree);
	tbr(childTree);
                                                             /* if the final tree cost is less than best cost computed */
	                                                     /* by SG algorithm, set this as best cost */
	if (childTree->costOfSubTask < bestCost) {
	  bestCost = childTree->costOfSubTask;
	}
	freeTree(childTree);

}
