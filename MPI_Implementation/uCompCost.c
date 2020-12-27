/**************************************************************************************************************************
 * DESCRIPTION: This function jumps from current node to its parent successively until it hits the root of the tree half. *
 * The parent holds information about its left and right children, these nodes are passed on to another function that     *
 * implements Fitch First Pass for measuring the cost.                                                                    *
 *                                                                                                                        *
 * INPUTS:                                                                                                                *
 * (1) Subtask containing the task to compute the length                                                                  *
 * (2) Index where node ID = -1 falls in the parent task of this subtask                                                  *
 * (3) Position where the new {internal node, taxa} pair was inserted in the parent subtask to get this child task        *
 * 	                                                                                                                  *
 * OUTPUTS:                                                                                                               *
 * (1) Returns cost of the tree. This value is also updated in the costOfThisSubTask element of the incoming subtask      *
 *************************************************************************************************************************/
#include "header.h"

void computeCost (struct subTask *sourceSubTask, int prevRhtSubTreeRootIndx, int pos) {

  int parentIndx, rChildIndx, lChildIndx, arrIndx, endIndx, strtIndx;
  struct subTaskArrayNode *leftChild, *rightChild, *parent;

                                                         /* Start from where the node was inserted in parent tree */
  strtIndx = pos;
                                                         /* End at the root of the subtree (if RHS: -1, if LHS: 0) */
  endIndx = (pos > prevRhtSubTreeRootIndx) ? sourceSubTask->rhtSubTreeRootIndx : 0; 

                                                         /* Loop until root of the subtree is encountered */
  arrIndx = strtIndx;
  while (arrIndx >= endIndx) {
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
                                                         /* If didnt reach the root of left/right subtree, jump to the */
                                                         /* parent of current node */
		if (arrIndx > 0 && arrIndx != sourceSubTask->rhtSubTreeRootIndx) {
		  parentIndx = sourceSubTask->node[arrIndx].parentIndx;
		  arrIndx = parentIndx;
		}
                                                         /* If reached the root of left/right subtree, break */
		else {
		  break;
		}
	}
  }
                                                         /* Both subtrees (left & right) are computed and cost kept */
                                                         /* in the roots of subtrees, compute cost of root node now */
  sourceSubTask->root.rChildIndx = sourceSubTask->rhtSubTreeRootIndx;
  sourceSubTask->root.lChildIndx = 0;
  performFitchOp(&sourceSubTask->node[0] , &sourceSubTask->root, &sourceSubTask->node[sourceSubTask->rhtSubTreeRootIndx]);
  sourceSubTask->costOfSubTask = sourceSubTask->root.costUntilThisNode;
}

/*************************************************************************************************************************
 * DESCRIPTION: This function implements fitch operation                                                                 *
 *                                                                                                                       *
 * INPUTS:                                                                                                               *
 * (1) Left child                                                                                                        *
 * (1) Right child                                                                                                       *
 * (1) Parent node                                                                                                       *
 * 	                                                                                                                 *
 * OUTPUTS:                                                                                                              *
 * (1) Computes cost of the node, and stores in the costUntilThisNode, costAtThisNode elements of the parent node        *
 *                                                                                                                       *
 *************************************************************************************************************************/
void performFitchOp(struct subTaskArrayNode *lChild, struct subTaskArrayNode *parent, struct subTaskArrayNode *rChild) {
   int site = 0;
   int length = 0;

#ifdef FITCH_OPERATION_VERBOSE
   /*
   printf("\n\t\t\t\t<><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
   printf("\t\t\t\t\t          FITCH OPERATION         \n");
   printf("\t\t\t\t<><><><><><><><><><><><><><><><><><><><><><><><><><>");
   printf("\n\t\t\t\t\tLeft child: %d, Parent: %d, Right Child: %d\n",lChild->id,parent->id,rChild->id);
   */
   printf("\t\t\t\t\t  Fitch Op {L:P:R}: {%4d:%4d:%4d}\n",lChild->id,parent->id,rChild->id);
#endif


   for (site=0; site < matrix.num_pars_inf_sites; site++) {
                                                               /* If parent needs a union of child states, incr length */
                                                               /* A intersection B is NULL, there is no common state */
      if ((lChild->stateVector[site] & rChild->stateVector[site])==0) {
        parent->stateVector[site] = lChild->stateVector[site] | rChild->stateVector[site];
	length = length + 1;
      }
                                                               /* Else, update parent's state to be a intersection */
                                                               /* A intersection B is not NULL, there is common state */
      else {
        parent->stateVector[site] = lChild->stateVector[site] & rChild->stateVector[site];
      }
   }
                                                               /* Update the cost of the MP tree until this node */
   parent->costAtThisNode = length;
   parent->costUntilThisNode = lChild->costUntilThisNode + rChild->costUntilThisNode + length;
}
