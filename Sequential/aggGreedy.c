/**************************************************************************************************************************
 * DESCRIPTION:                                                                                                           *
 *     This function sets the initial best cost. It uses smart greedy algorithm which finds the global best taxa, position*
 *     pair always, and adds it to the child tree                                                                         *
 *                                                                                                                        *
 * INPUTS:                                                                                                                *
 *     None                                                                                                               *
 *                                                                                                                        *
 * OUTPUTS:                                                                                                               *
 *     None                                                                                                               *
 *************************************************************************************************************************/

#include "header.h"

void getIstBestCost(void) {
	int taxaOver = 3;
	struct subTask *sourceTree, *childTree;
	struct maxMiniTaxaPosWLength mMTPWL;

#ifdef GET_IST_BESTCOST_VERBOSE
	printf("\n\t\t\t\t\t$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\t\t\t\t\t    IN GET IST BEST COST FUNCTION\n");
	printf("\t\t\t\t\t$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n");
#endif

	sourceTree = genIstSubTask();
#ifdef GET_IST_BESTCOST_VERBOSE
	printTree(sourceTree,-1);
#endif

	enqueueTaxa();
	mMTPWL = smartGreedy(sourceTree); 
	taxaOver += 1; 

	while (taxaOver < matrix.num_taxa) { 
#ifdef GET_IST_BESTCOST_VERBOSE
		printf("\n\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
		printf("\t\t\t\t\tADDING NEXT TAXA IN GET IST BEST COST\n");
		printf("\t\t\t\t\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
#endif
		childTree = genNextTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos); 
		freeTree(sourceTree);
		sourceTree = childTree;
#ifdef GET_IST_BESTCOST_VERBOSE
		printTree(sourceTree,mMTPWL.pos);
#endif
		mMTPWL = smartGreedy(sourceTree);
		taxaOver += 1; 
	}
	childTree = genNextTask(sourceTree, childTree, mMTPWL.taxa, mMTPWL.pos);
	freeTree(sourceTree);
	bestCost = mMTPWL.length;
	pushBestCostNode(childTree);
}

struct maxMiniTaxaPosWLength smartGreedy(struct subTask *sourceTree) { 
	struct maxMiniTaxaPosWLength prevmMTPWL, currmMTPWL; 
	struct taxaQ *tQNode;

	int taxaNum;

	prevmMTPWL.length = 1000000000;
	prevmMTPWL.taxa = 0;
	prevmMTPWL.pos = 0;
	currmMTPWL.length = 0; 
	currmMTPWL.taxa = 0; 
	currmMTPWL.pos = 0; 

#ifdef SMART_GREEDY_VERBOSE
	printf("\n\t*************************************************************************************************\n");
	printf("\n\t\t\t\t\t\tSMART GREEDY FUNCTION\n");
	printf("\n\t*************************************************************************************************\n\n");
#endif

	for (tQNode = taxaQHead; tQNode != NULL; tQNode=tQNode->next) { 
		taxaNum = tQNode->taxa; 
		currmMTPWL = getMinTree(sourceTree, taxaNum); 
#ifdef SMART_GREEDY_VERBOSE
		printf("\n\t\t->->->->->->RETAIN MINIMUM LENGTH OF ALL MINIMUM LENGTHS<-<-<-<-<-<-\n");
		printf("\t\tRecd {MinLen:%3d @ Pos:%3d w/ Taxa:%3d} vs. Prev {MinLen:%3d @ Pos:%3d w/ Taxa:%3d}\n\n",currmMTPWL.length, currmMTPWL.pos, taxaNum, prevmMTPWL.length, prevmMTPWL.pos, prevmMTPWL.taxa);
#endif
		if (currmMTPWL.length < prevmMTPWL.length) { 
#ifdef SMART_GREEDY_VERBOSE
			printf("\tCurr Min Lgt (%d) < Prev Min Lgt (%d)=> ",currmMTPWL.length,prevmMTPWL.length);
#endif
			prevmMTPWL.pos = currmMTPWL.pos; 
			prevmMTPWL.taxa = currmMTPWL.taxa; 
			prevmMTPWL.length = currmMTPWL.length; 
#ifdef SMART_GREEDY_VERBOSE
			printf("Updating: {CurrMinLen:%3d @ Pos:%3d w/ Taxon:%3d}\n\n",prevmMTPWL.length,prevmMTPWL.pos, prevmMTPWL.taxa);
#endif
		}
	}
	 
	deleteTaxa(prevmMTPWL.taxa); 

	return (prevmMTPWL); 
} 

