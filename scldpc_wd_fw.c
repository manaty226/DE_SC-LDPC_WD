#include<stdio.h>
#include<stdlib.h>
#include<math.h>

const int W = 16;		/* Window size */
const int dl = 3;		/* column weight */
const int dr = 6;		/* row weight */
const int GAMMA = 3;		/* Coupling width */
const double eps = 0.48812;	/* Channel erasure probability */
const int MAXITER = 100000;	/* Maximum number of decoding iterations */


/* Threshold analysis of SC-LDPC with windowed decoding */
int main()
{

  int i,j,k,l;
  double y[W + 2 * GAMMA];
  double y_res[W], chk, tmp;

  /* Initialisation */
  for(i = 0; i < W + 2 * GAMMA; i++) {
    if (i < GAMMA) {		
      y[i] = 0;
    }
    else {
      y[i] = 1;
    }
  }

  /* Decoding at the first window configuration */
  for(l = 0; l < MAXITER; l++) {
    /* Decoding at column nodes */
    for(i = 0; i < W; i++) {
      chk = 0;
      /* Decoding at row nodes */
      for(j = 0; j < GAMMA; j++) {
	tmp = 0;
	for(k = 0; k < GAMMA; k++) {
	  tmp += y[GAMMA + i + j - k];
	}
	chk += pow(1 - tmp / (double)GAMMA,dr-1);
      }
      y_res[i] = eps * pow(1 - chk/ (double)GAMMA, dl-1);
    }
    /* Update the message */
    for(i = 0; i < W; i++) {
      y[i + GAMMA] = y_res[i];
    }
  }
  
  /* Display erasure probability after decoding */
  for(i = 0; i < W + GAMMA * 2; i++) {
    printf("%d %.15f \n",i-2,y[i]);
  }

  return 0;
}
