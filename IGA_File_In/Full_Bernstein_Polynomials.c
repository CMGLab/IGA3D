#include <stdio.h>
#include <stdlib.h>
#include <math.h>

long Factorial(int n);

double Full_Berstein_Polynomials(double **B, int *p_1, int *p_2, int t, int elemN, int **IEN)
{
char temp[20];
int i,j,k,l,m,e;
int temp1;
double **B1 = malloc(elemN*sizeof(double*));
double **B2 = malloc(elemN*sizeof(double*));

int ii,nchoosek;
double tt,a,c,f;

for (e=0; e<elemN; e++){
	B1[e] = malloc((*(p_1+e)+1)*(t+1)*sizeof(double));
        for (j=0; j<t+1; j++){
		tt = ((j+1.0)-1.0)/t;
		for (k=0; k< (*(p_1+e)+1); k++){
			a = Factorial(*(p_1+e));
			c = Factorial(k);
			f = Factorial(*(p_1+e)-(k));
			nchoosek = a/(c*f);
			*(B1[e]+j*(*(p_1+e)+1)+(k)) = nchoosek*(pow(tt,k))*(pow((1-tt),(*(p_1+e)+1-(k+1))));
		}
	}
}

for (e=0; e<elemN; e++){
	B2[e] = malloc((*(p_2+e)+1)*(t+1)*sizeof(double));
        for (j=0; j<t+1; j++){
		tt = ((j+1.0)-1.0)/t;
		for (k=0; k< (*(p_2+e)+1); k++){
			a = Factorial(*(p_2+e));
			c = Factorial(k);
			f = Factorial(*(p_2+e)-(k));
			nchoosek = a/(c*f);
			*(B2[e]+j*(*(p_2+e)+1)+(k)) = nchoosek*(pow(tt,k))*(pow((1-tt),(*(p_2+e)+1-(k+1))));
		}
	}
}

for (e=0; e<elemN; e++){
	B[e] = malloc((*(p_1+e)+1)*(*(p_2+e)+1)*(t+1)*(t+1)*sizeof(double));
	for (i = 0; i<(*(p_1+e)+1); i++){
		for (j = 0; j<(*(p_2+e)+1); j++){
			for(k=0; k<t+1; k++){
				for (l=0; l<t+1; l++){
					ii = i*(*(p_2+e)+1)+j;
					*(B[e]+ii*(t+1)*(t+1)+k*(t+1)+l) = (*(B1[e]+k*(*(p_1+e)+1)+i))*(*(B2[e]+l*(*(p_2+e)+1)+j));
				}
			}
		}
	}
}

//free(B1);
//free(B2);
return **B;
}

