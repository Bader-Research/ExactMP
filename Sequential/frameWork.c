/*
 * DESCRIPTION:
 * This function sets the initial framework required to apply branch and bound
 * algorithm. The tasks it does are:
 * (1) Compute initial task with 3 taxa
 * (2) Preprocess the input data and assign State Encoding to the taxa
 * (3) Parsimonious Site Reordering
 * (4) Get initial best cost using SG algorithm
 * (5) Decide the taxa addition order using maxMini algorithm
 * (6) Push the first task to the task stack
 *                                     
 */

#include "header.h"

void setFrameWork(char *argv) {
  struct subTask *stRoot;
  struct taskNode *destNode;
  struct taxaQ *temp_node, *tQNode;

  int **rootTask=NULL,
      **childTask=NULL;

  numBestSolutions = 0;
  taskSTop = NULL;

  taxaQHead = NULL;
  taxaQTail = NULL;

  tAddOrdQHead = NULL;
  tAddOrdQTail = NULL;

  bestCostStackHead = NULL;

  taxaQueue = (int *)calloc(matrix.num_taxa, sizeof(int));
  assert(taxaQueue);

  printf("\nPreprocessing input...\n");
  preProcess(argv);
  printf("\nReordering sites...\n");
  reorderSites();

  //printf("\nApplying NJ Algorithm...\n");
  //nj();

  printf("\nApplying modified greedy algorithm...");
  fprintf(fp,"\nApplying modified greedy algorithm...");
  getIstBestCost();
  printf("\t\t[BEST SCORE %d]\n",bestCost+pNiCost);
  fprintf(fp,"\t\t[BEST SCORE %d]\n",bestCost+pNiCost);

  printf("\nDeciding taxa addition order...");
  fprintf(fp,"\nDeciding taxa addition order...");
  taxaAddOrder();
  
  printf("\nApplying eckDay greedy algorithm...");
  fprintf(fp,"\nApplying eckDay greedy algorithm...");
  eckDay();
  printf("\t\t[BEST SCORE %d]\n",bestCost+pNiCost);
  fprintf(fp,"\t\t[BEST SCORE %d]\n",bestCost+pNiCost);

  printf("\nApplying TBR on current best solution...");
  fprintf(fp,"\nApplying TBR on current best solution...");
  tbr(bestCostStackHead->bCNode);
  printf("\t[BEST SCORE %d]\n",bestCost+pNiCost);
  fprintf(fp,"\t[BEST SCORE %d]\n",bestCost+pNiCost);

  printf("\nApplying %d-level %d-arrangements/level randomized algorithm...\n",randOptLev,NUM_TBR_TREES);
  fprintf(fp,"\nApplying %d-level %d-arrangements/level randomized algorithm...\n",randOptLev,NUM_TBR_TREES);
  randOptimize();

  if (matrix.num_taxa > 7) {
      allocFitchArray();

      skipSize = 5 + matrix.num_pars_inf_sites,
      lhsOffSet = 4;
      rhsOffSet = 1+matrix.num_taxa;
      stack = (int *)calloc(matrix.num_taxa, sizeof(int));
      assert(stack);
      plbValues = (int *)calloc(matrix.num_taxa-6, sizeof(int));
      assert(plbValues);

      printf("\nComputing PLB for each level...\n");
      fprintf(fp,"\nComputing PLB for each level...\n");
      setPlbValues();
      printf("\nSetting up work space...\n");
      fprintf(fp,"\nSetting up work space...\n");
      allocMem4WorkSpace();
      getKthLevel(5);
  }
}

/*************************************************************************************************************************
 * DESCRIPTION:                                                                                                          *
 * This function prints the subtask in a human readable format. Nothing technical about it. No commenting required       *
 *************************************************************************************************************************/
void wr2File (struct subTask *sourceSubTask, FILE *fp, int bcNum) {
  int lc, rc, pr, i, arrIndx;
  struct subTask *sT;

  sT = sourceSubTask;

  fprintf(fp,"\n\t\t############################################################################################\n");
  fprintf(fp,"\t\t\t   Best Solution Number: %d {Tree Score (%d) = Task Cost (%d) + PNI Cost (%d)}\n",bcNum, sT->costOfSubTask+pNiCost, sT->costOfSubTask, pNiCost);
  fprintf(fp,"\t\t############################################################################################\n");
  fprintf(fp,"\n\tHighest Int Node: %d\t",sT->highestInternalNode);
  fprintf(fp,"Last taxon added: %d\t",sT->lastTaxaAdded);
  fprintf(fp,"Size of array : %d\t",sT->sizeOfSubTaskArray);
  fprintf(fp,"Cost of sub task: %d\n",sT->costOfSubTask);
  fprintf(fp,"\n\t\tARRAY ELEMENTS:\t");
  for (i=0; i < sT->sizeOfSubTaskArray; i++) {
	fprintf(fp,"[%d] ",sT->node[i].id);
  }
  fprintf(fp,"\n\n");

  fprintf(fp,"\t\t\tLCHILD\tNODE\tRCHILD\tPARENT\tCOST AT THIS NODE\t{STATE VECTOR}\n");
  for (arrIndx = sourceSubTask->sizeOfSubTaskArray-1; arrIndx >= 0; arrIndx = arrIndx-1) {
      rc = sT->node[arrIndx].rChildIndx;
      lc = sT->node[arrIndx].lChildIndx;
      pr = sT->node[arrIndx].parentIndx;
      if (sT->node[arrIndx].id > 0) {
        if (sT->node[arrIndx].parentIndx == ROOT) {
          fprintf(fp,"\t\t\t null\t%2d\t null\t root\t\t%2d\t\t",sT->node[arrIndx].id,sT->node[arrIndx].costAtThisNode);
	}
	else {
          fprintf(fp,"\t\t\t null\t%2d\t null\t%3d\t\t%2d\t\t",sT->node[arrIndx].id, sT->node[pr].id,sT->node[arrIndx].costAtThisNode);
	}
      }
      else {
        if (sT->node[arrIndx].parentIndx == ROOT) {
          fprintf(fp,"\t\t\t%3d\t%2d\t%3d\t root\t\t%2d\t\t",sT->node[lc].id, sT->node[arrIndx].id, sT->node[rc].id,sT->node[arrIndx].costAtThisNode);
	}
	else {
          fprintf(fp,"\t\t\t%3d\t%2d\t%3d\t%3d\t\t%2d\t\t",sT->node[lc].id, sT->node[arrIndx].id, sT->node[rc].id, sT->node[pr].id, sT->node[arrIndx].costAtThisNode);
	}
      }
      fprintf(fp,"{");
      for (i=0; i < matrix.num_pars_inf_sites; i++) {
         fprintf(fp,"%d ",sT->node[arrIndx].stateVector[i]);
      }
      fprintf(fp,"}\n");
  }
  fprintf(fp,"\t\t\t%3d\t%2d\t%3d\t%3d\t\t%2d\t\t",sT->node[sT->root.lChildIndx].id, sT->root.id, sT->node[sT->root.rChildIndx].id, sT->root.id, sT->root.costAtThisNode);
  fprintf(fp,"{");
  for (i=0; i < matrix.num_pars_inf_sites; i++) {
    fprintf(fp,"%d ",sT->root.stateVector[i]);
  }
  fprintf(fp,"}\n");

}
/*************************************************************************************************************************
 * DESCRIPTION:                                                                                                          *
 * This function prints the subtask in a human readable format. Nothing technical about it. No commenting required       *
 *************************************************************************************************************************/
void printTaxaStatus (struct taskNode *srcTaskNode) {
  struct taxaQ *tNode;
  //printf("\n\t\t\t\t\t---------------------------------------\n");
  //printf("\t\t\t\t\t    TAXA STATUS FOR ABOVE TREE \n");
  //printf("\t\t\t\t\t---------------------------------------\n");
  printf("\n\t\t\t\t    REMAINING TAXA FOR ABOVE TREE: ");
  for (tNode=srcTaskNode->taxaRemHead; tNode != NULL; tNode = tNode->next) {
    printf("%d ",tNode->taxa);
  }
  printf("\n\n\t\t\t\t    TAXA DONE FOR ABOVE TREE: ");
  for (tNode=srcTaskNode->taxaDoneHead; tNode != NULL; tNode = tNode->next) {
    printf("%d ",tNode->taxa);
  }
  printf("\n");
 // getchar();
}

/*************************************************************************************************************************
 * DESCRIPTION:                                                                                                          *
 * This function prints the subtask in a human readable format. Nothing technical about it. No commenting required       *
 *************************************************************************************************************************/
void printTree (struct subTask *sourceSubTask, int pos) {
  int lc, rc, pr, i, arrIndx;
  struct subTask *sT;
  
  sT = sourceSubTask;

  printf("\n\t\t\t\t\t---------------------------------------\n");
  printf("\t\t\t\t\t     TREE IN PREORDERED ARRAY FORM\n");
  printf("\t\t\t\t\t---------------------------------------\n");
  printf("\n\tHighest Int Node: %d\t",sT->highestInternalNode);
  printf("Last taxon added: %d\t",sT->lastTaxaAdded);
  if (pos != -1) {
    printf("in position: %d\t",pos);
  }
  printf("Size of array : %d\t",sT->sizeOfSubTaskArray);
  printf("Cost of sub task: %d\n",sT->costOfSubTask);
  //printf("\n\t\tARRAY ELEMENTS:\t");
  //for (i=0; i < sT->sizeOfSubTaskArray; i++) {
	//printf("[%d] ",sT->node[i].id);
  //}
  //printf("\n\n");

  printf("\t\t\tLCHILD\tNODE\tRCHILD\tPARENT\tCOST AT THIS NODE\t{STATE VECTOR}\n");
  for (arrIndx = sourceSubTask->sizeOfSubTaskArray-1; arrIndx >= 0; arrIndx = arrIndx-1) {
      rc = sT->node[arrIndx].rChildIndx;
      lc = sT->node[arrIndx].lChildIndx;
      pr = sT->node[arrIndx].parentIndx;
      if (sT->node[arrIndx].id > 0) {
        if (sT->node[arrIndx].parentIndx == ROOT) {
          printf("\t\t\t null\t%2d\t null\t root\t\t%2d\t\t",sT->node[arrIndx].id,sT->node[arrIndx].costAtThisNode);
	}
	else {
          printf("\t\t\t null\t%2d\t null\t%3d\t\t%2d\t\t",sT->node[arrIndx].id, sT->node[pr].id,sT->node[arrIndx].costAtThisNode);
	}
      }
      else {
        if (sT->node[arrIndx].parentIndx == ROOT) {
          printf("\t\t\t%3d\t%2d\t%3d\t root\t\t%2d\t\t",sT->node[lc].id, sT->node[arrIndx].id, sT->node[rc].id,sT->node[arrIndx].costAtThisNode);
	}
	else {
          printf("\t\t\t%3d\t%2d\t%3d\t%3d\t\t%2d\t\t",sT->node[lc].id, sT->node[arrIndx].id, sT->node[rc].id, sT->node[pr].id, sT->node[arrIndx].costAtThisNode);
	}
      }
      printf("{");
      for (i=0; i < matrix.num_pars_inf_sites; i++) {
         printf("%d ",sT->node[arrIndx].stateVector[i]);
      }
      printf("}\n");
  }
  printf("\t\t\t%3d\t%2d\t%3d\t%3d\t\t%2d\t\t",sT->node[sT->root.lChildIndx].id, sT->root.id, sT->node[sT->root.rChildIndx].id, sT->root.id, sT->root.costAtThisNode);
  printf("{");
  for (i=0; i < matrix.num_pars_inf_sites; i++) {
    printf("%d ",sT->root.stateVector[i]);
  }
  printf("}\n");
  //printf("\n\t\t---------------------------------------------------------------------------------\n");
  //printf("\t\t\t                            TREE OVER\n");
  //printf("\t\t---------------------------------------------------------------------------------\n");
}
