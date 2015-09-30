/*
 * IGA.cpp
 *
 *  Created on: Aug 19, 2015
 *      Author: Christopher
 */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <cblas.h>
#include <lapacke.h>
#include <petscksp.h>
#include "IGA.h"
using namespace std;

void Extract_Geometry(ASG_Element* mesh, int dim, int numElems, double** globalBezierPoints, double** globalBezierWeights) {

	// Loop through elements
	for (int i=0; i<numElems; i++) {

		int nloc = mesh[i].nloc;

		// Create array to hold control points and weights in projected space
		double projPts[nloc*(dim+1)];

		// Create array to hold Bezier control points and weights in projected space
		double bezPts[nloc*(dim+1)];

		// Multiply points by weights
		for (int j=0; j<nloc; j++) {
			for (int k=0; k<dim; k++)
				projPts[j*(dim+1)+k] = mesh[i].points[j*dim+k] * mesh[i].weights[j];
		}

		// Assign weights to final dimension of projected points array
		for (int j=0; j<nloc; j++)
			projPts[j*(dim+1)+dim] = mesh[i].weights[j];

		// Extract projected Bezier control points and weights
		cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,nloc,dim+1,nloc,1.0,mesh[i].extraction,nloc,projPts,dim+1,0.0,bezPts,dim+1);

		// Divide out weights and separate Bezier control points and weights:
		for (int j=0; j<nloc; j++) {
			for (int k=0; k<dim; k++)
				bezPts[j*(dim+1)+k] /= bezPts[j*(dim+1)+dim];
		}

		// Assign values to the arrays
		for (int j=0; j<nloc; j++) {
			for (int k=0; k<dim; k++)
				globalBezierPoints[i][j*dim+k] = bezPts[j*(dim+1)+k];
			globalBezierWeights[i][j] = bezPts[j*(dim+1)+dim];
		}

		mesh[i].bezierPoints = globalBezierPoints[i];
		mesh[i].bezierWeights = globalBezierWeights[i];
	}

	return;
}

void Quadrature_Data(int numQ, double* quadPts, double* quadWts) {

	switch (numQ) {
	case 1:
		quadPts[0] = 1./2.;
		quadWts[0] = 1.;
		break;
	case 2:
		quadPts[0] = 1./2.-1./2.*sqrt(1./3.);
		quadPts[1] = 1./2.+1./2.*sqrt(1./3.);
		quadWts[0] = 1./2.;
		quadWts[1] = 1./2.;
		break;
	case 3:
		quadPts[0] = 1./2.-1./2.*sqrt(3./5.);
		quadPts[1] = 1./2.;
		quadPts[2] = 1./2.+1./2.*sqrt(3./5.);
		quadWts[0] = 5./18.;
		quadWts[1] = 4./9.;
		quadWts[2] = 5./18.;
		break;
	case 4:
		quadPts[0] = 1./2.-1./2.*sqrt(3./7.-2./7.*sqrt(6./5.));
		quadPts[1] = 1./.2+1./2.*sqrt(3./7.-2./7.*sqrt(6./5.));
		quadPts[2] = 1./2.-1./2.*sqrt(3./7.+2./7.*sqrt(6./5.));
		quadPts[3] = 1./2.+1./2.*sqrt(3./7.+2./7.*sqrt(6./5.));
		quadWts[0] = 1./2.*(18.+sqrt(30.))/36.;
		quadWts[1] = 1./2.*(18.+sqrt(30.))/36.;
		quadWts[2] = 1./2.*(18.-sqrt(30.))/36.;
		quadWts[3] = 1./2.*(18.-sqrt(30.))/36.;
		break;
	}

	return;
}

void Bernstein_Basis_and_Deriv(double** ders, int p, double xi, int n) {

	// Adapted from Algorithm A2.3 in The NURBS Book

	double ndu[p+1][p+1];
	ndu[0][0] = 1;

	double left[p+1];
	double right[p+1];
	double d, j1, j2, saved, temp;
	int b, pk, rk, s1, s2;
	double a[2][2];

	for (int j=1; j<=p; j++) {
		left[j] = xi - 0;
		right[j] = 1- xi;
		saved = 0;

		for (int r=0; r<j; r++) {
			ndu[j][r] = right[r+1] + left[j-r];
			temp = ndu[r][j-1]/ndu[j][r];

			ndu[r][j] = saved+right[r+1]*temp;
			saved = left[j-r]*temp;
		}

		ndu[j][j] = saved;
	}

	for (int j=0; j<=p; j++)
		ders[0][j] = ndu[j][p];

	for (int r=0; r<=p; r++) {
		s1 = 0;
		s2 = 1;
		a[0][0] = 1;

		for (int k=1; k<=n; k++) {
			d = 0.;
			rk = r-k;
			pk = p-k;

			if (r >=k ) {
				a[s2][0] = a[s1][0]/ndu[pk+1][rk];
				d = a[s2][0]*ndu[rk][pk];
			}

			if (rk >= -1)
				j1 = 1;
			else
				j1 = -rk;

			if (r-1 <= pk)
				j2 = k-1;
			else
				j2 = p-r;

			for (int j=j1; j<=j2; j++) {
				a[s2][j] = (a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
				d += a[s2][j]*ndu[rk+j][pk];
			}

			if (r <= pk) {
				a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
				d += a[s2][k]*ndu[r][pk];
			}

			ders[k][r] = d;
			b = s1;
			s1 = s2;
			s2 = b;
		}
	}

	double r = (double) p;
	for (int k=1; k<=n; k++) {
		for (int j=0; j<=p; j++) {
			ders[k][j] *= r;
		}
		r *= (p-k);
	}

	return;
}

void Bernstein_Basis_and_Derivs(double** B, int p1, int p2, double xi_1, double xi_2, int n) {

	double** N;
	N = new double* [n+1];
	for (int i=0; i<n+1; i++)
		N[i] = new double [p1+1];

	double** M;
	M = new double* [n+1];
	for (int i=0; i<n+1; i++)
		M[i] = new double [p2+1];

	Bernstein_Basis_and_Deriv(N,p1,xi_1,n);
	Bernstein_Basis_and_Deriv(M,p2,xi_2,n);

	int a = 0;

	// Indexing gets weird here due to changes in convention
	for (int i=0; i<p1+1; i++) {
		for (int j=0; j<p2+1; j++) {
			//Basis functions
			B[a][0] = N[0][i]*M[0][j];
			//First derivatives
			B[a][1] = N[1][i]*M[0][j];
			B[a][2] = N[0][i]*M[1][j];
			if (n==2) {
				//Second derivatives
				B[a][3] = N[0][i]*M[2][j];
				B[a][4] = N[1][i]*M[1][j];
				B[a][5] = N[2][i]*M[0][j];
			}
			a++;
		}
	}

	// Free variables
	for (int i=0; i<n+1; i++)
		delete[] N[i];
	delete[] N;

	for (int i=0; i<n+1; i++)
		delete[] M[i];
	delete[] M;

	return;
}

void Full_Bernstein_Basis_and_Derivs(double*** fullBasis, int nloc, int p, int n_q, int n, double* quadPts) {

	int index;

	// Define total number of functions to store (including all derivatives)
	int numFuncts = ((n+1)*(n+2))/2;

	// Create variable to hold Bernstein basis and derivatives at a given point
	double** B;
	B = new double* [nloc];
	for (int i=0; i<nloc; i++)
		B[i] = new double [numFuncts];

	// Loop through quadrature points and calculate Bernstein basis and derivatives at each point
	for (int i=0; i<n_q; i++) {
		for (int j=0; j<n_q; j++) {
			Bernstein_Basis_and_Derivs(B,p,p,quadPts[i],quadPts[j],n);
			index = i*n_q+j;

			for (int k=0; k<nloc; k++) {
				for (int l=0; l<numFuncts; l++)
					fullBasis[k][l][index] = B[k][l];
			}
		}
	}

	// Free variables
	for (int i=0; i<nloc; i++)
		delete[]  B[i];
	delete[] B;

	return;
}

void Boundary_Bernstein_Basis_and_Derivs(double*** fullBasis, int nloc, int p, int n_q, int n, double* quadPts) {


	int index;

	// Define total number of functions to store (including all derivatives)
	int numFuncts = ((n+1)*(n+2))/2;

	// Create variable to hold Bernstein basis and derivatives at a given point
	double** B;
	B = new double* [nloc];
	for (int i=0; i<nloc; i++)
		B[i] = new double [numFuncts];

	// Loop through boundaries and calculate Bernstein basis and derivatives at each point along each boundary
	double xi_t1, xi_t2;

	for (int i=0; i<4; i++) {
		for (int j=0; j<n_q; j++) {

			switch(i) {
			case 0:
				xi_t1 = quadPts[j];
				xi_t2 = 0;
				break;
			case 1:
				xi_t1 = 1;
				xi_t2 = quadPts[j];
				break;
			case 2:
				xi_t1 = quadPts[j];
				xi_t2 = 1;
				break;
			case 3:
				xi_t1 = 0;
				xi_t2 = quadPts[j];
				break;
			}

			Bernstein_Basis_and_Derivs(B,p,p,xi_t1,xi_t2,n);
			index = i*n_q+j;

			for (int k=0; k<nloc; k++) {
				for (int l=0; l<numFuncts; l++)
					fullBasis[k][l][index] = B[k][l];
			}
		}
	}

	// Free variables
	for (int i=0; i<nloc; i++)
		delete[] B[i];
	delete[] B;

	return;
}

void ShapeFunction(ASG_Element element, double basis[], double R[], double dRdx[], double d2Rdx2[], double x[], double J[], int flag) {

	int n_loc = element.nloc;

	for (int i=0; i<n_loc; i++)
		R[i] = 0;
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<2; j++)
			//dRdx[i][j] = 0;
			dRdx[i*2+j] = 0;
	}

	x[0] = 0;
	x[1] = 0;

	double wb = 0;
	double dwbdxi[2] = {0, 0};
	double dRdxi[n_loc][2];
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<2; j++)
			dRdxi[i][j] = 0;
	}

	double dxdxi[2][2];
	for (int i=0; i<2; i++) {
		for (int j=0; j<2; j++)
			dxdxi[i][j] = 0;
	}

	double dxidx[2][2];
	for (int i=0; i<2; i++) {
		for (int j=0; j<2; j++)
			dxidx[i][j] = 0;
	}

	double B[n_loc];
	for (int i=0; i<n_loc; i++)
		B[i] = basis[i*6];

	double dBdxi[n_loc][2];
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<2; j++)
			dBdxi[i][j] = basis[i*6+(j+1)];
	}

	double d2Bdxi2[n_loc][3];
	if (flag == 1) {
		for (int i=0; i<n_loc; i++) {
			for (int j=0; j<3; j++)
				d2Bdxi2[i][j] = basis[i*6+(j+3)];
		}
	}

	double d2wbdxi2[3];
	for (int i=0; i<3; i++)
		d2wbdxi2[i] = 0;

	double d2Rdxi2[n_loc][3];
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<3; j++)
			d2Rdxi2[i][j] = 0;
	}

	double d2xdxi2[2][3];
	for (int i=0; i<2; i++) {
		for (int j=0; j<3; j++)
			d2xdxi2[i][j] = 0;
	}

	double bp[n_loc][2];

	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<2; j++)
			bp[i][j] = element.bezierPoints[i*2+j];
	}

	//weighting and derivatives
	for (int a=0; a<n_loc; a++){
		wb += element.bezierWeights[a]*B[a];
		dwbdxi[0] += element.bezierWeights[a]*dBdxi[a][0];
		dwbdxi[1] += element.bezierWeights[a]*dBdxi[a][1];

		if (flag == 1) {
			d2wbdxi2[0] += element.bezierWeights[a]*d2Bdxi2[a][0];
			d2wbdxi2[1] += element.bezierWeights[a]*d2Bdxi2[a][1];
			d2wbdxi2[2] += element.bezierWeights[a]*d2Bdxi2[a][2];
			}
	}

	//basis functions and parametric derivatives
	double C_operators[n_loc][n_loc];

	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<n_loc; j++)
			C_operators[i][j] = element.extraction[i*n_loc+j];
	}

	for (int a=0; a<n_loc; a++) {
		for (int b=0; b<n_loc; b++) {
			R[a] += element.weights[a] * C_operators[a][b]*B[b]/wb;
			for (int i=0; i<2; i++)
				dRdxi[a][i] += element.weights[a] * C_operators[a][b] * (dBdxi[b][i]/wb - dwbdxi[i]*B[b]/(wb*wb));

			if (flag == 1) {
				d2Rdxi2[a][0] += element.weights[a]*C_operators[a][b]*((d2Bdxi2[b][0]/wb - dBdxi[b][0]*dwbdxi[0]/(wb*wb)) - (dwbdxi[0]*(dBdxi[b][0]/(wb*wb) - 2*dwbdxi[0]*B[b]/(wb*wb*wb)) + d2wbdxi2[0]*B[b]/(wb*wb)));
				d2Rdxi2[a][1] += element.weights[a]*C_operators[a][b]*((d2Bdxi2[b][1]/wb - dBdxi[b][0]*dwbdxi[1]/(wb*wb)) - (dwbdxi[0]*(dBdxi[b][1]/(wb*wb) - 2*dwbdxi[1]*B[b]/(wb*wb*wb)) + d2wbdxi2[1]*B[b]/(wb*wb)));
				d2Rdxi2[a][2] += element.weights[a]*C_operators[a][b]*((d2Bdxi2[b][2]/wb - dBdxi[b][1]*dwbdxi[1]/(wb*wb)) - (dwbdxi[1]*(dBdxi[b][1]/(wb*wb) - 2*dwbdxi[1]*B[b]/(wb*wb*wb)) + d2wbdxi2[2]*B[b]/(wb*wb)));
			}
		}
	}

	//physical space quantities
	for (int a=0; a<n_loc; a++) {
		for (int i=0; i<2; i++)
			x[i] += element.bezierWeights[a]*bp[a][i]*B[a]/wb;
		for (int i=0; i<2; i++) {
			for (int j=0; j<2; j++)
				dxdxi[i][j] += element.bezierWeights[a]*bp[a][i]*(dBdxi[a][j]/wb - dwbdxi[j]*B[a]/(wb*wb));
		}

		if (flag == 1) {
			for (int i=0; i<2; i++) {
				d2xdxi2[i][0] += element.bezierWeights[a]*bp[a][i]*((d2Bdxi2[a][0]/wb - dBdxi[a][0]*dwbdxi[0]/(wb*wb)) - (dwbdxi[0]*(dBdxi[a][0]/(wb*wb) - 2*dwbdxi[0]*B[a]/(wb*wb*wb)) + d2wbdxi2[0]*B[a]/(wb*wb)));
				d2xdxi2[i][1] += element.bezierWeights[a]*bp[a][i]*((d2Bdxi2[a][1]/wb - dBdxi[a][0]*dwbdxi[1]/(wb*wb)) - (dwbdxi[0]*(dBdxi[a][1]/(wb*wb) - 2*dwbdxi[1]*B[a]/(wb*wb*wb)) + d2wbdxi2[1]*B[a]/(wb*wb)));
				d2xdxi2[i][2] += element.bezierWeights[a]*bp[a][i]*((d2Bdxi2[a][2]/wb - dBdxi[a][1]*dwbdxi[1]/(wb*wb)) - (dwbdxi[1]*(dBdxi[a][1]/(wb*wb) - 2*dwbdxi[1]*B[a]/(wb*wb*wb)) + d2wbdxi2[2]*B[a]/(wb*wb)));
			}
		}
	}

	for (int i=0; i<2; i++) {
		for (int j=0; j<2; j++)
			J[i*2+j] = dxdxi[i][j];
	}

	double detJ = J[0]*J[3] - J[1]*J[2];
	dxidx[0][0] = J[3]/detJ;
	dxidx[0][1] = -J[1]/detJ;
	dxidx[1][0] = -J[2]/detJ;
	dxidx[1][1] = J[0]/detJ;

	//physical derivatives
	for (int a=0; a<n_loc; a++) {
		for (int i=0; i<2; i++) {
			for (int j=0; j<2; j++)
				//dRdx[a][i] += dRdxi[a][j]*dxidx[j][i];
				dRdx[a*2+i] += dRdxi[a][j]*dxidx[j][i];
		}
	}

	if (flag == 1) {
		double A[9];

		A[0] = dxdxi[0][0]*dxdxi[0][0];
		A[1] = 2*dxdxi[0][0]*dxdxi[1][0];
		A[2] = dxdxi[1][0]*dxdxi[1][0];
		A[3] = dxdxi[0][0]*dxdxi[0][1];
		A[4] = dxdxi[0][0]*dxdxi[1][1] + dxdxi[0][1]*dxdxi[1][0];
		A[5] = dxdxi[1][0]*dxdxi[1][1];
		A[6] = dxdxi[0][1]*dxdxi[0][1];
		A[7] = 2*dxdxi[0][1]*dxdxi[1][1];
		A[8] = dxdxi[1][1]*dxdxi[1][1];

		double b[3*n_loc];
		double LHS[3*n_loc];

		for (int i=0; i<n_loc; i++) {
			//b[i*3] = dRdx[i][0]*d2xdxi2[0][0] + dRdx[i][1]*d2xdxi2[1][0];
			//b[i*3+1] = dRdx[i][0]*d2xdxi2[0][1] + dRdx[i][1]*d2xdxi2[1][1];
			//b[i*3+2] = dRdx[i][0]*d2xdxi2[0][2] + dRdx[i][1]*d2xdxi2[1][2];
			b[i*3] = dRdx[i*2]*d2xdxi2[0][0] + dRdx[i*2+1]*d2xdxi2[1][0];
			b[i*3+1] = dRdx[i*2]*d2xdxi2[0][1] + dRdx[i*2+1]*d2xdxi2[1][1];
			b[i*3+2] = dRdx[i*2]*d2xdxi2[0][2] + dRdx[i*2+1]*d2xdxi2[1][2];

			LHS[i*3] = d2Rdxi2[i][0] - b[i*3];
			LHS[i*3+1] = d2Rdxi2[i][1] - b[i*3+1];
			LHS[i*3+2] = d2Rdxi2[i][2] - b[i*3+2];
		}

		int ipiv[3];

		LAPACKE_dgetrf(LAPACK_ROW_MAJOR,3,3,A,3,ipiv);
		LAPACKE_dgesv(LAPACK_ROW_MAJOR,3,n_loc,A,3,ipiv,LHS,n_loc);

		for (int i=0; i<n_loc; i++) {
			for (int j=0; j<3; j++)
				//d2Rdx2[i][j] = LHS[i*3+j];
				d2Rdx2[i*3+j] = LHS[i*3+j];
		}
	}

	return;
}



double Grid_Spacing(double bezierPoints[], int p) {

	//TODO Generalize to 3D

	int nloc = (p+1)*(p+1);
	double h_grid, xRange, yRange;
	int dim = 2;

	xRange = bezierPoints[nloc*dim-2] - bezierPoints[0];
	yRange = bezierPoints[nloc*dim-1] - bezierPoints[1];

	if (xRange <= yRange)
		h_grid = xRange;
	else
		h_grid = yRange;

	return h_grid;
}



void Element_Assembly(ASG_Element element, int* BC, double* g, double* k, double* f, Mat K, Vec F) {

	int nloc = element.nloc;
	int i, j;

	for (int a=0; a<nloc; a++){
		i = element.connectivity[a];
		if (BC[i] == 0) {
			for (int b=0; b<nloc; b++) {
				j = element.connectivity[b];
				if (BC[j] == 1) {
					f[a] = f[a] - k[a*nloc+b]*g[j];
					k[a*nloc+b] = 0;
				}
			}
		}
		else {
			for (int b=0; b<nloc; b++) {
				j = element.connectivity[b];
				k[a*nloc+b] = 0;
			}
			f[a] = 0;
		}
	}

	MatSetValues(K,element.nloc,element.connectivity,element.nloc,element.connectivity,k,ADD_VALUES);
	VecSetValues(F,element.nloc,element.connectivity,f,ADD_VALUES);

	/*for (int a=0; a<nloc; a++){
		i = element.connectivity[a];
		if (BC[i] == 0) {
			for (int b=0; b<nloc; b++) {
				j = element.connectivity[b];
				if (BC[j] == 0) {
					//MatSetValue(K,i,j,k[a*nloc+b],ADD_VALUES);
					Kvals[i*numNodes+j] += k[a*nloc+b];
				} else {
					//VecSetValue(F,i,-k[a*nloc+b]*g[j],ADD_VALUES);
					Fvals[i] -= k[a*nloc+b]*g[j];
				}
			}
			//VecSetValue(F,i,f[a],ADD_VALUES);
			Fvals[i] += f[a];
		}
	}*/

	return;
}

double Get_Local_Value(double vi[], double basis[], int nloc) {

	double locVal = 0;

	for (int i=0; i<nloc; i++)
		locVal += vi[i]*basis[i];

	return locVal;
}

void Extract_Concentration(ASG_Element* mesh, double*** fullBasis, int numElems, double* plotPts, int numPlotPts, double* d, double* X, double* Y, double* U) {

	// TODO: Generalize to handle curves in 3D

	// Assume all elements have same nloc
	int nloc = mesh[0].nloc;
	double xi1, xi2, u;

	double basis[nloc*6];

	double R[nloc];

	//double** dRdx;
	//dRdx = new double* [nloc];
	//for (int j=0; j<nloc; j++)
	//	dRdx[j] = new double [2];
	double dRdx[nloc*2];

	//double** d2Rdx2;
	//d2Rdx2 = new double* [nloc];
	//for (int j=0; j<nloc; j++)
	//	d2Rdx2[j] = new double [3];
	double d2Rdx2[nloc*3];

	double x[2];

	double J[2*2];

	int flag = 0;

	for (int e=0; e<numElems; e++) {
		for (int i1=0; i1<numPlotPts; i1++) {
			xi1 = plotPts[i1];

			for (int i2=0; i2<numPlotPts; i2++) {
				xi2 = plotPts[i2];

				for (int j=0; j<nloc; j++) {
					for (int k=0; k<3; k++)
						basis[j*6+k] = fullBasis[j][k][i1*numPlotPts+i2];
				}

				ShapeFunction(mesh[e],basis,R,dRdx,d2Rdx2,x,J,flag);

				u = 0;
				for (int a=0; a<nloc; a++)
					u += R[a]*d[mesh[e].connectivity[a]];

				// Uncommenting the lines below causes the correct (x,y) values to be written
				/*cout << "yo 1854" << endl;
				for (int i=0; i<nloc; i++)
					cout << R[i] << ", ";
				cout << endl;*/

				X[e*numPlotPts*numPlotPts+i1*numPlotPts+i2] = x[0];
				Y[e*numPlotPts*numPlotPts+i1*numPlotPts+i2] = x[1];
				U[e*numPlotPts*numPlotPts+i1*numPlotPts+i2] = u;
			}
		}
	}

	// Free variables
	//for (int j=0; j<nloc; j++)
	//	delete[] dRdx[j];
	//delete[] dRdx;

	//for (int j=0; j<nloc; j++)
	//	delete[] d2Rdx2[j];
	//delete[] d2Rdx2;

	return;
}

void Subelement_Connectivity(int* subConnectivity, int numPlotPts, int numSubElems) {

	// TODO: Adjust this to incorporate subdivision of subelements
	int e;

	for (int i=0; i<(numPlotPts-1); i++) {
		for (int j=0; j<(numPlotPts-1); j++) {
			e = i*(numPlotPts-1)+j;
			subConnectivity[e] = i*numPlotPts+j;
			subConnectivity[e+numSubElems] = i*numPlotPts+j+1;
			subConnectivity[e+2*numSubElems] = i*numPlotPts+j+numPlotPts;
			subConnectivity[e+3*numSubElems] = i*numPlotPts+j+numPlotPts+1;
		}
	}

	return;
}
