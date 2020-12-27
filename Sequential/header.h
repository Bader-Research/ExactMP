#include <math.h>
#include <time.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

#define ROOT -9999999
#define LEAF -8888888
#define TAXA_NAME_LENGTH 100
#define MAXMINI
//#define NUM_TBR_TREES 600
#define THREE_OPTIMIZATION

#define TRUE         1
#define FALSE        0
#define MAX_RAND     2147483647
#define DBL_INF      (double)MAX_RAND
#define MAXSTATES    32
#define ALLTREES
//#define ONETREE

//#define MAXMAX
//#define TBR_VERY_VERBOSE
//#define TBR_VERIFY
//#define TBR_VERBOSE
//#define BnB_VERBOSE
//#define REORDER_VERBOSE
//#define MAX_MINI_VERBOSE
//#define BnB_VERY_VERBOSE
//#define BnB_VERY_VERY_VERBOSE
#define GET_ALL_SOLUTIONS   
//#define GET_ONE_BEST_SOLUTION
//#define PREPROCESS_VERBOSE
//#define SMART_GREEDY_VERBOSE
//#define GET_MIN_TREE_VERBOSE
//#define TAXA_ADD_ORDER_VERBOSE
//#define FITCH_OPERATION_VERBOSE
//#define GET_IST_BESTCOST_VERBOSE

#define keepTrees 100

FILE *fp;

int bestCost,
    initBestCost,
    pNiCost,
    *tbrBest,
    numPNISites,
    numConstSites,
    numBestSolutions,
    MPScoreFromNJ,
    NJUniqStates,                        /* unique states in sequence matrix */
    *stack,
    skipSize,
    lhsOffSet,
    rhsOffSet,
    ****workSpace,
    ***solQ,
    solQCtr,
    *plbValues,
    *taxaQueue,
    randOptFlag,
    randOptLev,
    NUM_TBR_TREES,
    ***fitchArray,
    **wosac;                             // Work Space counter, a 2D arrray
                                         // where 1st dimension is for the row
                                         // in work space. 1st element in 2nd 
                                         // dimension is ctr in that row and 
                                         // 2nd element in that dimension is 
                                         // size of that row

struct bestCostStack {
  struct subTask *bCNode;
  struct bestCostStack *down;
};

struct bestCostStack *bestCostStackHead;

struct integerStack {
   int id;
   struct integerStack *down;
};

struct integerStack *stackTop;

struct {
   char *uniq_chars;                           /* Characters representing the unique states */
   double **taxon_table;                       /* Index r,c holds evolutionary distance between ith & jth taxa */
} NJ;

struct {
   int  num_taxa;                              /* Total number of input taxa (rows) */
   int  num_sites;                             /* Total number of input sites (cols) */
   char **taxons;                              /* Index 0 of this array holds the name of 0th taxon */
   char **taxon_table;                         /* Index r,c of this array holds the rth taxa's cth character */
   int  **state_encoding;                      /* Index [r][c] of this array holds encoding of taxon_table[r][c] */
   int  **reord_sites_enc;                     /* Reordered array containing only parsimony informative sites */
   int  num_pars_inf_sites;                    /* Total number of parsimony informative sites */
   int  *num_st_at_site;                       /* Index 0 of this holds the # unique states of site 0 */
   char **unique_states;                       /* Holds unique states in a site */
   int  **unique_states_reps;                  /* Index 0 of this array holds # reps of index 0 of unique_states array */
} matrix;

struct subTaskArrayNode {                      // (6 + matrix.num_pars_inf_sites) integers
  int id;                                      // 1 integer
  int parentIndx;                              // 1 integer
  int lChildIndx;                              // 1 integer
  int rChildIndx;                              // 1 integer
  int costAtThisNode;                          // 1 integer
  int costUntilThisNode;                       // 1 integer
  int *stateVector;                            // matrix.num_pars_inf_sites integers
};

struct subTask {                               // (sizeOfSubTaskArray+1) X (6+matrix.num_pars_inf_sites) + 4
  int sizeOfSubTaskArray;                      // 1 integer
  int lastTaxaAdded;                           // 1 integer
  int highestInternalNode;                     // 1 integer
  int costOfSubTask;                           // 1 integer
  int rhtSubTreeRootIndx;                      // 0 integer                        NO LONGER EXISTANT IN 1D ARRAY
  struct subTaskArrayNode root;                // 6 + matrix.num_pars_inf_sites
  struct subTaskArrayNode *node;               // sizeOfSubTaskArray X (6 + matrix.num_pars_inf_sites)
};

struct taskNode {
  struct taxaQ *taxaRemHead;
  struct taxaQ *taxaRemTail;
  struct taxaQ *taxaDoneHead;
  struct taxaQ *taxaDoneTail;
  struct subTask *partialTask;
};

struct taskStack {
  struct taskNode *tNode;
  struct taskStack *down;
};

struct taxaQ {
  int taxa;
  struct taxaQ *next;
};

struct taskStack *taskSTop;

struct taxaQ *taxaQHead;
struct taxaQ *taxaQTail;

struct taxaQ *tAddOrdQHead;
struct taxaQ *tAddOrdQTail;

struct maxMiniTaxaPosWLength {
   int pos;
   int taxa;
   int length;
};

void   checkUsage(int, char **);
void   tbr(struct subTask *);
int    pop(void);
void   push(int);
void   showBits(int);
void   deleteTaxa(int );
void   freeMemory(void );
void   preProcess(char *);
void   reorderSites(void );
void   taxaAddOrder(void);
void   getIstBestCost(void);
void   branchAndBound(void);
struct taskNode *popTask(void);
void   flushBestCostStack(void);
int    maskStateInPos(int , int );
struct subTask* dequeueTask(void);
void   freeTree(struct subTask *);
void   pushTask(struct taskNode *);
struct subTask *genIstSubTask(void);
void   setFrameWork(char *);
void randOptimize(void);
void   enqueueTaxa(void);
void   enqueueTask(struct subTask *);
int    dequeueTaxaQ(struct taskNode *);
void   printTree(struct subTask *, int);
void   freeStackNode(struct taskNode *);
int    purdomsLowerBnd(struct taskNode *);
void   pushBestCostNode(struct subTask *);
void   subTaskPrint(struct subTask *, int);
void   printTaxaStatus (struct taskNode *);
void   wr2File(struct subTask *, FILE *, int);
struct taskNode *copyTaxaQ(struct taskNode *);
void   computeCost(struct subTask *, int, int);
void   setParentChildIndices (struct subTask *);
void   printBestSolutions(int, int, int, int, FILE *);
struct maxMiniTaxaPosWLength maxMini(struct subTask *);
struct subTask *genNextTask(struct subTask *, struct subTask *, int , int );
struct maxMiniTaxaPosWLength smartGreedy(struct subTask *);
struct maxMiniTaxaPosWLength getMinTree(struct subTask *, int );
struct maxMiniTaxaPosWLength getMaxTree(struct subTask *, int );
void performFitchOp(struct subTaskArrayNode *, struct subTaskArrayNode *, struct subTaskArrayNode *);
struct subTask *genNextEckDayTask(struct subTask *, struct subTask *, int , int );
void eckDay(void);
void computeEckDayCost (struct subTask *);
struct maxMiniTaxaPosWLength getEckDayMinTree(struct subTask *, int );
void TBReckDay(int **, int );
void searchOptimize(void);
void correctEckDayIndices (struct subTask *);
void TBRmaxMin(void);
void nQ4RandOpt(int **, int );
void TBRSmartGreedy(void);
void nj(void);
void neighborJoin(void);
void clearMemory(void);
void enqueueAddOrdTaxa(int ); 
int **genIstOptSubTask(int **);
int numNodesInRow(int );
void getKthLevel(int );
void allocMem4WorkSpace(void);
void clearWorkSpace(void);
void setPlbValues(void);
void bnb(void);
void genTaskLeft(int **, int **, int, int);
void genTaskRight(int **, int **, int, int);
void freeRandTaxaQ(void);
void displayTree(int **, int , int );
void displayTreeInNexusFormat(int **, int , int );
