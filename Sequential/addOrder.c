/*************************************************************************************************************************
 * DESCRIPTION: This function generates the taxa addition order using max mini algorithm.                                *
 *                                                                                                                       *
 * INPUTS: None                                                                                                          *
 *                                                                                                                       *
 * OUTPUTS: Queue containing the addition order of taxa, writes directly into global queue. Doesnt return anything       *
 ************************************************************************************************************************/
#include "header.h"

void taxaAddOrder(void) {
	int taxaOver = 3;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;
	struct taxaQ *t;

	sourceTree = genIstSubTask();
                                                             /* enqueue all taxa to be checked for addition order */
                                                             /* this will be from 4 to NUMBER of TAXA */
	enqueueTaxa();
                                                             /* apply maxmini, & enqueue first addition order taxa */
	mMTPWL = maxMini(sourceTree); 
	enqueueAddOrdTaxa(mMTPWL.taxa); 
	taxaOver += 1; 
                                                             /* for all taxa remaining to be enqueued in addition order */
	while (taxaOver < matrix.num_taxa) { 
#ifdef TAXA_ADD_ORDER_VERBOSE
		printf("\t\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~\n");
		printf("\t\t\t\t\t\tADDING NEXT TAXA ORDER\n");
		printf("\t\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~");
#endif
                                                             /* generate next task with previous maxmini taxa, position */
		childTree = genNextTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos); 
                                                             /* free parent task, it is not needed anymode */
		freeTree(sourceTree);
		sourceTree = childTree;
#ifdef TAXA_ADD_ORDER_VERBOSE
		printTree(sourceTree,mMTPWL.pos);
#endif
                                                             /* apply maxmini, & enqueue next addition order taxa */
		mMTPWL = maxMini(sourceTree);
		enqueueAddOrdTaxa(mMTPWL.taxa); 
		taxaOver += 1; 
	}
                                                             /* if the final tree cost is less than best cost computed */
	                                                     /* by SG algorithm, set this as best cost */
	if (mMTPWL.length < bestCost) {
	  bestCost = mMTPWL.length;
	  flushBestCostStack();
  	  childTree = genNextTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos);
	  freeTree(sourceTree);
	  pushBestCostNode(childTree);
	}
	else {
	  freeTree(sourceTree);
	}

	printf("\t\t\t[BEST SCORE %d]\n",bestCost+pNiCost);
	fprintf(fp,"\t\t\t[BEST SCORE %d]\n",bestCost+pNiCost);

        printf("\tTaxa addition order: 1 2 3 ");
        fprintf(fp,"\tTaxa addition order: 1 2 3 ");
	taxaQueue[0] = 1;
	taxaQueue[1] = 2;
	taxaQueue[2] = 3;
        for (taxaOver=3, t=tAddOrdQHead; t != NULL; t = t->next, taxaOver=taxaOver+1) {
           taxaQueue[taxaOver] = t->taxa;
           printf("%d ",t->taxa);
           fprintf(fp,"%d ",t->taxa);
        }
        printf("\n");
        fprintf(fp,"\n");

}

/*************************************************************************************************************************
 * DESCRIPTION: This function implements maxmini algorithm                                                               *
 *                                                                                                                       *
 * INPUTS: Sub task on which to apply maxmini algorithm                                                                  *
 *                                                                                                                       *
 * OUTPUTS: structure containing the taxon to be added next, position where it is to be added and tree lenght            *
 ************************************************************************************************************************/
struct maxMiniTaxaPosWLength maxMini(struct subTask *sourceTree) { 
	struct maxMiniTaxaPosWLength prevmMTPWL, currmMTPWL; 
	struct taxaQ *tQNode;

	int taxaNum;

	prevmMTPWL.length = 0;
	prevmMTPWL.taxa = 0;
	prevmMTPWL.pos = 0;
	currmMTPWL.length = 0; 
	currmMTPWL.taxa = 0; 
	currmMTPWL.pos = 0; 

#ifdef MAX_MINI_VERBOSE
	printf("\n\t\t*************************************************************************************************\n");
	printf("\n\t\t\t\t\t\t\tMAXMINI FUNCTION\n");
	printf("\n\t\t*************************************************************************************************\n");
#endif
                                                             /* for all taxa remaining to be added */
	for (tQNode = taxaQHead; tQNode != NULL; tQNode=tQNode->next) { 
		taxaNum = tQNode->taxa; 
                                                             /* get tree with minimum length with this taxa added to the */
                                                             /* subtask under view (check for all positions) */
#ifdef MAXMINI
		currmMTPWL = getMinTree(sourceTree, taxaNum); 
#endif
#ifdef MAXMAX
		currmMTPWL = getMaxTree(sourceTree, taxaNum); 
#endif
#ifdef MAX_MINI_VERBOSE
		printf("\n\t\t->->->->->->RETAIN MAXIMUM LENGTH OF ALL MINIMUM LENGTHS<-<-<-<-<-<-\n");
                printf("\t\tCurr {MinLen:%3d @ Pos:%3d w/ Taxa:%3d} vs. Prev {MinLen:%3d @ Pos:%3d w/ Taxa:%3d}\n",currmMTPWL.length, currmMTPWL.pos, currmMTPWL.taxa, prevmMTPWL.length, prevmMTPWL.pos, prevmMTPWL.taxa);
#endif
                                                             /* this length > previous length, this is maximum, keep it */
		if (currmMTPWL.length > prevmMTPWL.length) { 
#ifdef MAX_MINI_VERBOSE
			printf("\t\tCurr Min Lgt (%d) > Prev Min Lgt (%d)=> ",currmMTPWL.length,prevmMTPWL.length);
#endif
			prevmMTPWL.pos = currmMTPWL.pos; 
			prevmMTPWL.taxa = currmMTPWL.taxa; 
			prevmMTPWL.length = currmMTPWL.length; 
#ifdef MAX_MINI_VERBOSE
                        printf("Update Curr Max Len to {MaxLen:%3d @ Pos:%3d w/ Taxa:%3d}\n", prevmMTPWL.length, prevmMTPWL.pos, prevmMTPWL.taxa);
#endif
		}
	}
                                                             /* remove taxa added to taxa addition order */
	deleteTaxa(prevmMTPWL.taxa);
	 
	return (prevmMTPWL); 
} 
 
/**************************************************************************************************************************
 * DESCRIPTION: This function returns the tree with minimum length, with the position where a given taxa should be added  *
 *              This checks the given taxa at all positions and returns the minimum tree. If there is a tie, it just picks*
 *              up the first one and ignores the rest                                                                     *
 *                                                                                                                        *
 * INPUTS: Sub task on which to apply maxmini algorithm, and taxa which is to be added                                    *
 *                                                                                                                        *
 * OUTPUTS: Structure containing the taxon to be added next, position where it is to be added and tree length             *
 *************************************************************************************************************************/
struct maxMiniTaxaPosWLength getMinTree(struct subTask *sourceTree, int taxaNum) { 
	struct subTask *childTree; 
	int currentLength = 0; 
	int prevLength = 0; 
	 
	struct maxMiniTaxaPosWLength mMTPWL; 
	 
	int pos;
	int currMinPos; 
	int currMinTaxa; 
	int arraySize = sourceTree->sizeOfSubTaskArray; 
	 
#ifdef GET_MIN_TREE_VERBOSE
	printf("\n\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	printf("\t\t\t\t CHECKING MINIMUM TREE LENGTH FOR ALL POSITIONS FOR TAXA: %d\n",taxaNum);
	printf("\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
#endif
                                                             /* generate first case, adding taxa in position 0 */
	childTree = genNextTask(sourceTree, childTree, taxaNum, 0); 
#ifdef GET_MIN_TREE_VERBOSE
	printTree(childTree, 0);
#endif
                                                             /* update initial tree length to the length in position 0 */
	prevLength = childTree->costOfSubTask;
	freeTree(childTree);

	currMinPos = 0;
	currMinTaxa = taxaNum; 

                                                             /* for all remaining positions */
	for (pos=1; pos < arraySize; pos++) { 
                                                             /* do not try position before internal node -1 */
                                                             /* this is done to ensure that the integrity of tree */
		if (sourceTree->node[pos].id != -1) { 
			childTree = genNextTask(sourceTree,childTree, taxaNum, pos); 
#ifdef GET_MIN_TREE_VERBOSE
	                printTree(childTree, pos);
#endif
			currentLength = childTree->costOfSubTask; 
#ifdef GET_MIN_TREE_VERBOSE
			printf("\n\t\t->->->->->->RETAIN MINIMUM LENGTH FROM CURRENT & PREVIOUS LENGTHS<-<-<-<-<-<-\n");
                        printf("\n\t\t\tCurr {MinLen:%3d @ Pos:%3d w/ Taxa:%3d} vs. Prev {MinLen:%3d @ Pos:%3d w/ Taxa:%3d}\n\n",currentLength, pos, taxaNum, prevLength, currMinPos, currMinTaxa);
#endif
                                                             /* if this length is less than previous length, keep it */
			if (currentLength < prevLength) {
#ifdef GET_MIN_TREE_VERBOSE
				printf("\t\tCurr Lgt (%d) < Prev Lgt (%d)=> ",currentLength,prevLength);
#endif
				prevLength = currentLength; 
				currMinPos = pos; 
				currMinTaxa = taxaNum; 
#ifdef GET_MIN_TREE_VERBOSE
				printf("Update Curr Lgt to: %d & keeping Taxa: %d, Pos: %d\n\n",currentLength,currMinTaxa,currMinPos);
#endif
			} 

			freeTree(childTree);
		}
	}
		 
                                                             /* return position & tree with minimum length */
	mMTPWL.taxa = currMinTaxa;
	mMTPWL.pos = currMinPos;
	mMTPWL.length = prevLength;
		 
	return (mMTPWL); 
}

/**************************************************************************************************************************
 * DESCRIPTION: This function returns the tree with maximum length, with the position where a given taxa should be added  *
 *              This checks the given taxa at all positions and returns the maximum tree. If there is a tie, it just picks*
 *              up the first one and ignores the rest                                                                     *
 *                                                                                                                        *
 * INPUTS: Sub task on which to apply maxmax algorithm, and taxa which is to be added                                     *
 *                                                                                                                        *
 * OUTPUTS: Structure containing the taxon to be added next, position where it is to be added and tree length             *
 *************************************************************************************************************************/
struct maxMiniTaxaPosWLength getMaxTree(struct subTask *sourceTree, int taxaNum) { 
	struct subTask *childTree; 
	int currentLength = 0; 
	int prevLength = 0; 
	 
	struct maxMiniTaxaPosWLength mMTPWL; 
	 
	int pos;
	int currMinPos; 
	int currMinTaxa; 
	int arraySize = sourceTree->sizeOfSubTaskArray; 
	 
#ifdef GET_MIN_TREE_VERBOSE
	printf("\n\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	printf("\t\t\t\t CHECKING MINIMUM TREE LENGTH FOR ALL POSITIONS FOR TAXA: %d\n",taxaNum);
	printf("\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
#endif
                                                             /* generate first case, adding taxa in position 0 */
	childTree = genNextTask(sourceTree, childTree, taxaNum, 0); 
#ifdef GET_MIN_TREE_VERBOSE
	printTree(childTree, 0);
#endif
                                                             /* update initial tree length to the length in position 0 */
	prevLength = childTree->costOfSubTask;
	freeTree(childTree);

	currMinPos = 0;
	currMinTaxa = taxaNum; 

                                                             /* for all remaining positions */
	for (pos=1; pos < arraySize; pos++) { 
                                                             /* do not try position before internal node -1 */
                                                             /* this is done to ensure that the integrity of tree */
		if (sourceTree->node[pos].id != -1) { 
			childTree = genNextTask(sourceTree,childTree, taxaNum, pos); 
#ifdef GET_MIN_TREE_VERBOSE
	                printTree(childTree, pos);
#endif
			currentLength = childTree->costOfSubTask; 
#ifdef GET_MIN_TREE_VERBOSE
			printf("\n\t\t->->->->->->RETAIN MAXIMUM LENGTH FROM CURRENT & PREVIOUS LENGTHS<-<-<-<-<-<-\n");
                        printf("\n\t\t\tCurr {MinLen:%3d @ Pos:%3d w/ Taxa:%3d} vs. Prev {MinLen:%3d @ Pos:%3d w/ Taxa:%3d}\n\n",currentLength, pos, taxaNum, prevLength, currMinPos, currMinTaxa);
#endif
                                                             /* if this length is more than previous length, keep it */
			if (currentLength > prevLength) { 
#ifdef GET_MIN_TREE_VERBOSE
				printf("\t\tCurr Lgt (%d) > Prev Lgt (%d)=> ",currentLength,prevLength);
#endif
				prevLength = currentLength; 
				currMinPos = pos; 
				currMinTaxa = taxaNum; 
#ifdef GET_MIN_TREE_VERBOSE
				printf("Update Curr Lgt to: %d & keeping Taxa: %d, Pos: %d\n\n",currentLength,currMinTaxa,currMinPos);
#endif
			} 

			freeTree(childTree);
		}
	}
		 
                                                             /* return position & tree with minimum length */
	mMTPWL.taxa = currMinTaxa;
	mMTPWL.pos = currMinPos;
	mMTPWL.length = prevLength;
		 
	return (mMTPWL); 
}
