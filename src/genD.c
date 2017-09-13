//
//  genD.c
//  
//
//  Created by Marco Enea on 02/10/15.
//
//

#include "genD.h"
#include <math.h>
#include <Rmath.h>
#include <R.h>

void genD(int *ncols, int *d){
    int i, j;
    int k=0;
    int nc = *ncols;
    for (i = 1; i < *ncols; i++)
    {
      for (j = 0; j <= (i-1); j++)
      {
        d[k+i] = 1;
        d[k+j] = -1;
        k = k + nc;
      }
    } 
}
