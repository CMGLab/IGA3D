/*
 * VMS.cpp
 *
 *  Created on: Aug 19, 2015
 *      Author: Christopher
 */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <petscksp.h>
#include <cblas.h>
#include <lapacke.h>
#include <petscksp.h>
#include "VMS.h"
#include "IGA.h"
using namespace std;

void ShapeFine(ASG_Element element, int pf, double** basis, double** dRfdx, double* J) {

	int n_loc = (pf+1)*(pf+1);

	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<2; j++)
			dRfdx[i][j] = 0;
	}

	double dxidx[2][2];

	double dBdxi[n_loc][2];
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<2; j++)
			dBdxi[i][j] = basis[i][j+1];
	}

	double detJ = J[0]*J[3] - J[1]*J[2];
	dxidx[0][0] = J[3]/detJ;
	dxidx[0][1] = -J[1]/detJ;
	dxidx[1][0] = -J[2]/detJ;
	dxidx[1][1] = J[0]/detJ;

	// physical derivatives
	for (int a=0; a<n_loc; a++) {
		for (int i=0; i<2; i++) {
			for (int j=0; j<2; j++)
				dRfdx[a][i] += dBdxi[a][j]*dxidx[j][i];
		}
	}
}

void Solve_Fine_Scales(ASG_Element element, BC_Struct elementBC, Prob_Params probParams, double*** fullBasis, double*** fullBasisFine, double*** fullBoundaryBasis, double*** fullBoundaryBasisFine, int n_q, double* quadWts, int pf, double* d, double* d_f) {

	int n_loc = element.nloc;
	int n_locf = (pf+1)*(pf+1);

	double k_fc[n_locf*n_loc];
	for (int i=0; i<n_locf; i++) {
		for (int j=0; j<n_loc; j++)
			k_fc[i*n_loc+j] = 0;
	}

	double k_ff[n_locf*n_locf];
	for (int i=0; i<n_locf; i++) {
		for (int j=0; j<n_locf; j++)
			k_ff[i*n_locf+j] = 0;
	}

	double f_f[n_locf];
	for (int i=0; i<n_locf; i++)
		f_f[i] = 0;

	double phi_c[n_loc];

	double basis[n_loc*6];

	double** basisFine;
	basisFine = new double* [n_locf];
	for (int j=0; j<n_locf; j++)
		basisFine[j] = new double [3];

	double boundaryBasis[n_loc*6];

	double** boundaryBasisFine;
	boundaryBasisFine = new double* [n_locf];
	for (int j=0; j<n_locf; j++)
		boundaryBasisFine[j] = new double [3];

	double R[n_loc];

	//double** dRdx;
	//dRdx = new double* [n_loc];
	//for (int j=0; j<n_loc; j++)
	//	dRdx[j] = new double [2];
	double dRdx[n_loc*2];

	//double** d2Rdx2;
	//d2Rdx2 = new double* [n_loc];
	//for (int j=0; j<n_loc; j++)
	//	d2Rdx2[j] = new double [3];
	double d2Rdx2[n_loc*3];

	double x[2];

	double J[2*2];

	int flag = 1;

	double Rf[n_locf];

	double** dRfdx;
	dRfdx = new double* [n_locf];
	for (int j=0; j<n_locf; j++)
		dRfdx[j] = new double [2];

	double f_loc, Jdet, k_art, k_loc;
	double v_loc[2];
	double u_1[n_loc];
	double u_2[n_loc];

	double advGalIntegral, advTerm, artDiffIntegral, diffGalIntegral, diffTerm;

	// Calculate grid spacing
	double h_grid;
	h_grid = Grid_Spacing(element.bezierPoints,element.p);

	// Loop through quadrature points

	for (int q1=0; q1<n_q; q1++) {
		for (int q2=0; q2<n_q; q2++) {

			for (int j=0; j<n_loc; j++) {
				for (int k=0; k<6; k++)
					basis[j*6+k] = fullBasis[j][k][q1*n_q+q2];
			}

			for (int j=0; j<n_locf; j++) {
				for (int k=0; k<3; k++)
					basisFine[j][k] = fullBasisFine[j][k][q1*n_q+q2];
			}

			ShapeFunction(element,basis,R,dRdx,d2Rdx2,x,J,flag);
			ShapeFine(element,pf,basisFine,dRfdx,J);
			for (int i=0; i<n_locf; i++)
				Rf[i] = basisFine[i][0];

			Jdet = J[0]*J[3] - J[1]*J[2];

			k_loc = Get_Local_Value(probParams.kappa,R,n_loc);
			k_art = Get_Local_Value(probParams.kappaArt,R,n_loc);
			f_loc = Get_Local_Value(probParams.f,R,n_loc);

			for (int i=0; i<n_loc; i++) {
				u_1[i] = probParams.u[i];
				u_2[i] = probParams.u[n_loc+i];
			}

			v_loc[0] = Get_Local_Value(u_1,R,n_loc);
			v_loc[1] = Get_Local_Value(u_2,R,n_loc);

			for (int af=0; af<n_locf; af++) {
				for (int b=0; b<n_loc; b++){
					//advTerm = (v_loc[0]*dRdx[b][0] + v_loc[1]*dRdx[b][1])*Rf[af]*quadWts[q1]*quadWts[q2]*Jdet;
					//diffTerm = -k_loc*Rf[af]*(d2Rdx2[b][0] + d2Rdx2[b][2])*quadWts[q1]*quadWts[q2]*Jdet;
					advTerm = (v_loc[0]*dRdx[b*2] + v_loc[1]*dRdx[b*2+1])*Rf[af]*quadWts[q1]*quadWts[q2]*Jdet;
					diffTerm = -k_loc*Rf[af]*(d2Rdx2[b*3] + d2Rdx2[b*3+2])*quadWts[q1]*quadWts[q2]*Jdet;
					k_fc[af*n_loc+b] += advTerm + diffTerm;

					phi_c[b] = d[element.connectivity[b]];
				}
			}

			for (int af=0; af<n_locf; af++) {
				for (int bf=0; bf<n_locf; bf++) {
					diffGalIntegral = k_loc*(dRfdx[af][0]*dRfdx[bf][0] + dRfdx[af][1]*dRfdx[bf][1])*quadWts[q1]*quadWts[q2]*Jdet;
					advGalIntegral = -Rf[bf]*(v_loc[0]*dRfdx[af][0] + v_loc[1]*dRfdx[af][1])*quadWts[q1]*quadWts[q2]*Jdet;
					artDiffIntegral = k_art*(dRfdx[af][0]*dRfdx[bf][0] + dRfdx[af][1]*dRfdx[bf][1])*quadWts[q1]*quadWts[q2]*Jdet;
					k_ff[af*n_locf+bf] += diffGalIntegral + advGalIntegral + artDiffIntegral;
				}

				f_f[af] += Rf[af]*f_loc*quadWts[q1]*quadWts[q2]*Jdet;
			}
		}
	}

	//Boundary conditions
	double t_t[2];
	double tangent[2];
	double t_n[2];
	double C, h_loc, tangentNorm, j_boundary;
	double advNeumBC, diffNeumBC, diffTerm1, diffTerm2, penaltyTerm;

	for (int side=0; side<4; side++) {
		for (int q=0; q<n_q; q++) {
			switch(side) {
			case 0:
				t_t[0] = -1;
				t_t[1] = 0;
				break;
			case 1:
				t_t[0] = 0;
				t_t[1] = -1;
				break;
			case 2:
				t_t[0] = 1;
				t_t[1] = 0;
				break;
			case 3:
				t_t[0] = 0;
				t_t[1] = 1;
				break;
			}

			for (int j=0; j<n_loc; j++) {
				for (int k=0; k<3; k++)
					boundaryBasis[j*6+k] = fullBoundaryBasis[j][k][side*n_q+q];
			}

			for (int j=0; j<n_locf; j++) {
				for (int k=0; k<3; k++)
					boundaryBasisFine[j][k] = fullBoundaryBasisFine[j][k][side*n_q+q];
			}

			flag = 0;

			ShapeFunction(element,boundaryBasis,R,dRdx,d2Rdx2,x,J,flag);
			ShapeFine(element,pf,boundaryBasisFine,dRfdx,J);
			for (int i=0; i<n_locf; i++)
				Rf[i] = boundaryBasisFine[i][0];

			tangent[0] = J[0]*t_t[0] + J[1]*t_t[1];
			tangent[1] = J[2]*t_t[0] + J[3]*t_t[1];

			tangentNorm = sqrt(tangent[0]*tangent[0] + tangent[1]*tangent[1]);

			t_n[0] = -tangent[1]/tangentNorm;
			t_n[1] = tangent[0]/tangentNorm;

			j_boundary = tangentNorm;

			k_loc = Get_Local_Value(probParams.kappa,R,n_loc);

			for (int i=0; i<n_loc; i++) {
				u_1[i] = probParams.u[i];
				u_2[i] = probParams.u[n_loc+i];
			}

			v_loc[0] = Get_Local_Value(u_1,R,n_loc);
			v_loc[1] = Get_Local_Value(u_2,R,n_loc);

			C = probParams.C;

			for (int af=0; af<n_locf; af++) {
				for (int bf=0; bf<n_locf; bf++) {

					if (v_loc[0]*t_n[0] + v_loc[1]*t_n[1] > 0)
						advTerm = Rf[af]*Rf[bf]*(v_loc[0]*t_n[0] + v_loc[1]*t_n[1])*quadWts[q]*j_boundary;
					else
						advTerm = 0;

					diffTerm1 = -k_loc*(t_n[0]*dRfdx[bf][0] + t_n[1]*dRfdx[bf][1])*Rf[af]*quadWts[q]*j_boundary;
					diffTerm2 = -k_loc*(t_n[0]*dRfdx[af][0] + t_n[1]*dRfdx[af][1])*Rf[bf]*quadWts[q]*j_boundary;
					penaltyTerm = k_loc*C/h_grid*Rf[af]*Rf[bf]*quadWts[q]*j_boundary;

					k_ff[af*n_locf+bf] += advTerm + diffTerm1 + diffTerm2 + penaltyTerm;
				}
			}
		}
	}

	// Invert k_ff
	int ipiv[n_locf];
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR,n_locf,n_locf,k_ff,n_locf,ipiv);
	LAPACKE_dgetri(LAPACK_ROW_MAJOR,n_locf,k_ff,n_locf,ipiv);

	// Calculate d_f

	cblas_dgemv(CblasRowMajor,CblasNoTrans,n_locf,n_loc,-1.0,k_fc,n_loc,phi_c,1,1.0,f_f,1);
	cblas_dgemv(CblasRowMajor,CblasNoTrans,n_locf,n_locf,1.0,k_ff,n_locf,f_f,1,0.0,d_f,1);

	// Free variables
	for (int j=0; j<n_locf; j++)
		delete[] basisFine[j];
	delete[] basisFine;

	for (int j=0; j<n_locf; j++)
		delete[] boundaryBasisFine[j];
	delete[] boundaryBasisFine;

	//for (int j=0; j<n_loc; j++)
	//	delete[] dRdx[j];
	//delete[] dRdx;

	//for (int j=0; j<n_loc; j++)
	//	delete[] d2Rdx2[j];
	//delete[] d2Rdx2;

	for (int j=0; j<n_locf; j++)
		delete[] dRfdx[j];
	delete[] dRfdx;

	return;
}

void Extract_Fine_Scales(double*** fullBasisFine, int numElems, int numPlotPts, int pf, double* d, double* U) {

	// TODO: Generalize to handle curves in 3D

	// Assume all elements have same nlocf
	int nlocf = (pf+1)*(pf+1);
	double u;

	double Rf[nlocf];

	double basisFine[nlocf][3];

	for (int e=0; e<numElems; e++) {
		for (int i1=0; i1<numPlotPts; i1++) {
			for (int i2=0; i2<numPlotPts; i2++) {

				for (int j=0; j<nlocf; j++)
					Rf[j] = fullBasisFine[j][0][i1*numPlotPts+i2];

				u = 0;
				for (int a=0; a<nlocf; a++)
					u += Rf[a]*d[e*nlocf+a];

				U[e*numPlotPts*numPlotPts+i1*numPlotPts+i2] = u;
			}
		}
	}

	return;
}
