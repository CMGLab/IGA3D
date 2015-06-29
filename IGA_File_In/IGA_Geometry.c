#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double Full_Berstein_Polynomials(double **B, int *p_1, int *p_2, int t, int elemN, int **IEN);


int IGA_Geometry(int *nloc, double **Pb, double **wb, double **C, int *p_1, int *p_2, int **IEN, double **B, int *n , int d)
{
char temp[20];
int i,j,k,l,m,e;
int temp1;
int nodeN,elemN,t=1;
/*-------------READ in input data for geometry____________*/
d=3;
FILE *ifp;
ifp = fopen("Membrane.iga","r");
fscanf(ifp,"%s %s",temp,temp);
fscanf(ifp,"%s %d",temp,&nodeN);
fscanf(ifp,"%s %d",temp,&elemN);
double *P = malloc(2*nodeN*d*sizeof(double));
double *w = malloc(nodeN*sizeof(double));

for (i=0; i<nodeN; i++){
        fscanf(ifp,"%s",temp);
        for (j=0; j<d; j++){
                fscanf(ifp,"%lf",(P+i*d+j));
	}
	fscanf(ifp,"%lf\n",(w+i));       
}

for (i=0; i<elemN; i++){
	fscanf(ifp,"%s %d %d %d\n",temp,n+i,p_1+i,p_2+i);
	*(nloc+i)=(*(p_1+i)+1)*(*(p_2+i)+1);
        IEN[i] = malloc((*(nloc+i))*sizeof(int));
	C[i] = malloc((*(nloc+i))*(*(n+i))*sizeof(double));
	for (j=0; j<*(n+i); j++){
		fscanf(ifp,"%d",(IEN[i]+j));
	}
        fscanf(ifp,"\n");
	for (k=0; k< *(n+i) ; k++){
		for (l=0; l<nloc[i]  ; l++){
			fscanf(ifp,"%lf",(C[i]+k*(*(nloc+i))+l));
		}
	fscanf(ifp,"\n");
	} 	
}
fclose(ifp);


/*__________________Perform Extraction____________________*/

double *Pp = malloc(nodeN*(d+1)*sizeof(double));
for (i=0; i<nodeN; i++){
	for (j=0; j<d; j++){
		*(Pp+i*(d+1)+j) = (*(P+i*d+j))*(*(w+i));
	}
}

for (i=0; i<nodeN; i++){
	*(Pp+i*(d+1)+d) = *(w+i);
}

double **pbp = malloc(elemN*sizeof(double*));
for (i=0; i < elemN; i ++){
pbp[i] = malloc((*(nloc+i))*(d+1)*sizeof(double));
	for ( j=0; j< (d+1);j++){
		for (k=0; k<nloc[i]; k++){
			for (l =0; l<n[i]; l++){
				*(pbp[i]+j*nloc[i]+k) = *(pbp[i]+j*nloc[i]+k)  +  (*(C[i]+l*nloc[i]+k))  *  (*(Pp+(*(IEN[i]+l))*(d+1)+j));
			}
		}
	}
}

for (i=0; i<elemN; i++){
	wb[i] = malloc(nloc[i]*sizeof(double));
	for(k=0; k<nloc[i]; k++){
		*(wb[i]+k) = *(pbp[i]+d*nloc[i]+k);
	}
}

for (i=0; i<elemN; i++){
	Pb[i] = malloc(nloc[i]*d*sizeof(double));
	for (j=0; j<d; j++){
		for( k=0; k< nloc[i]; k++){
			*(Pb[i]+j*nloc[i]+k) = *(pbp[i]+j*nloc[i]+k)/(*(wb[i]+k));
		}
	}
} 

/*___________Compute Element Basis Functions______________*/

Full_Berstein_Polynomials(B,p_1,p_2,t,elemN,IEN);
double ww;

/*___________Compute Surface Points for Geometry__________*/

double *S = malloc(elemN*d*(t+1)*(t+1)*sizeof(double));
for (j=0; j<d; j++){
	for (e=0; e<elemN; e++){
		for (k=0; k<t+1; k++){
			for (i=0; i<t+1; i++){
				
				ww=0;
				for(l=0; l<nloc[e]; l++){
					ww = ww + (*(wb[e]+l)) * (*(B[e]+l*(t+1)*(t+1)+k*(t+1)+i));
				}
									
				for(m=0; m<nloc[e]; m++){
					*(S+j*elemN*(t+1)*(t+1)+e*(t+1)*(t+1)+k*(t+1)+i) =  (*(S+j*elemN*(t+1)*(t+1)+e*(t+1)*(t+1)+k*(t+1)+i)) +
					(((*(Pb[e]+j*nloc[e]+m))  *    (*(wb[e]+m))   *  (*(B[e]+m*(t+1)*(t+1)+k*(t+1)+i)))  / ww);
			
				}
			}			
		}
	}
}

/*___________Output Surface to be Plotted_________________*/

FILE *ofp;
ofp = fopen("surface.txt","w");

for (e=0; e<elemN; e++){
	for (i=0; i<t+1; i++){
		for(j=0; j<t+1; j++){
			fprintf(ofp,"%lf \t %lf \t %lf\n",*(S+0*elemN*(t+1)*(t+1)+e*(t+1)*(t+1)+i*(t+1)+j),*(S+1*elemN*(t+1)*(t+1)+e*(t+1)*(t+1)+i*(t+1)+j),
			*(S+2*elemN*(t+1)*(t+1)+e*(t+1)*(t+1)+i*(t+1)+j));
		}
	}
}

fclose(ofp);

FILE *oofp;
oofp = fopen("param.txt","w");
fprintf(oofp,"%d \t %d\n",t,elemN);

//free(S);
//free(P);//
//free(w);
//free(Pp);
//for(i=0; i<elemN; i++){
//	free(pbp[i]);
//}
//free(pbp);
return 0;
}

