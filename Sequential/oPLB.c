#include "header.h"

void setPlbValues(void) {
  int done, rem, statesOver, statesRem, site, plbLen, plbIndx, mask, 
      bpos, maskRem, maskOver;
  

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
  /*
  for (plbIndx=7; plbIndx < matrix.num_taxa; plbIndx++) {
	  printf("PLB{%d}: %d\n",plbIndx-7, plbValues[plbIndx-7]);
  }
  getchar();
  */
}
