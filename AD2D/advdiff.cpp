/*
 * advdiff.cpp
 *
 *  Created on: Aug 19, 2015
 *      Author: Christopher
 */

#include <cmath>
#include <cblas.h>
#include <lapacke.h>
#include "advdiff.h"
#include "IGA.h"
#include "VMS.h"
using namespace std;

void Element_Formation(ASG_Element element, BC_Struct elementBC, Prob_Params probParams, double*** fullBasis, double*** fullBasisFine, double*** fullBoundaryBasis, double*** fullBoundaryBasisFine, int n_q, double* quadWts, int pf, double* k_e, double* f_e) {

	int n_loc = element.nloc;
	int n_locf = (pf+1)*(pf+1);

	double k_cc[n_loc*n_loc];
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<n_loc; j++)
			k_cc[i*n_loc+j] = 0;
	}

	double f_c[n_loc];
	for (int i=0; i<n_loc; i++)
		f_c[i] = 0;

	double k_fc[n_locf*n_loc];
	for (int i=0; i<n_locf; i++) {
		for (int j=0; j<n_loc; j++)
			k_fc[i*n_loc+j] = 0;
	}

	double k_cf[n_loc*n_locf];
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<n_locf; j++)
			k_cf[i*n_locf+j] = 0;
	}

	double k_ff[n_locf*n_locf];
	for (int i=0; i<n_locf; i++) {
		for (int j=0; j<n_locf; j++)
			k_ff[i*n_locf+j] = 0;
	}

	double f_f[n_locf];
	for (int i=0; i<n_locf; i++)
		f_f[i] = 0;

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

			for (int a=0; a<n_loc; a++) {
				for (int b=0; b<n_loc; b++) {
					//diffGalIntegral = k_loc*(dRdx[a][0]*dRdx[b][0] + dRdx[a][1]*dRdx[b][1])*quadWts[q1]*quadWts[q2]*Jdet;
					//advGalIntegral = -R[b]*(v_loc[0]*dRdx[a][0] + v_loc[1]*dRdx[a][1])*quadWts[q1]*quadWts[q2]*Jdet;
					diffGalIntegral = k_loc*(dRdx[a*2]*dRdx[b*2] + dRdx[a*2+1]*dRdx[b*2+1])*quadWts[q1]*quadWts[q2]*Jdet;
					advGalIntegral = -R[b]*(v_loc[0]*dRdx[a*2] + v_loc[1]*dRdx[a*2+1])*quadWts[q1]*quadWts[q2]*Jdet;
					k_cc[a*n_loc+b] += diffGalIntegral + advGalIntegral;
				}
				diffTerm = R[a]*f_loc*quadWts[q1]*quadWts[q2]*Jdet;
				f_c[a] += diffTerm;
			}

			for (int af=0; af<n_locf; af++) {
				for (int b=0; b<n_loc; b++){
					//advTerm = (v_loc[0]*dRdx[b][0] + v_loc[1]*dRdx[b][1])*Rf[af]*quadWts[q1]*quadWts[q2]*Jdet;
					//diffTerm = -k_loc*Rf[af]*(d2Rdx2[b][0] + d2Rdx2[b][2])*quadWts[q1]*quadWts[q2]*Jdet;
					advTerm = (v_loc[0]*dRdx[b*2] + v_loc[1]*dRdx[b*2+1])*Rf[af]*quadWts[q1]*quadWts[q2]*Jdet;
					diffTerm = -k_loc*Rf[af]*(d2Rdx2[b*3] + d2Rdx2[b*3+2])*quadWts[q1]*quadWts[q2]*Jdet;
					k_fc[af*n_loc+b] += advTerm + diffTerm;
				}
			}

			for (int a=0; a<n_loc; a++) {
				for (int bf=0; bf<n_locf; bf++) {
					//k_cf[a*n_locf+bf] += -(v_loc[0]*dRdx[a][0] + v_loc[1]*dRdx[a][1])*Rf[bf]*quadWts[q1]*quadWts[q2]*Jdet;
					k_cf[a*n_locf+bf] += -(v_loc[0]*dRdx[a*2] + v_loc[1]*dRdx[a*2+1])*Rf[bf]*quadWts[q1]*quadWts[q2]*Jdet;
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

			if (elementBC.Neumann[side] == 1) {

				h_loc = elementBC.h[side];

				for (int a=0; a<n_loc; a++) {
					for (int b=0; b<n_loc; b++) {

						if (v_loc[0]*t_n[0] + v_loc[1]*t_n[1] > 0)
							advNeumBC = R[a]*R[b]*(v_loc[0]*t_n[0] + v_loc[1]*t_n[1])*quadWts[q]*j_boundary;
						else
							advNeumBC = 0;

						k_cc[a*n_loc+b] += advNeumBC;
					}

					diffNeumBC = R[a]*h_loc*quadWts[q]*j_boundary;
					f_c[a] += diffNeumBC;
				}
			}
		}
	}

	// Invert k_ff
	int ipiv[n_locf];
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR,n_locf,n_locf,k_ff,n_locf,ipiv);
	LAPACKE_dgetri(LAPACK_ROW_MAJOR,n_locf,k_ff,n_locf,ipiv);

	// Calculate k_e and f_e

	double temp[n_loc*n_locf];

	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n_loc,n_locf,n_locf,1.0,k_cf,n_locf,k_ff,n_locf,0.0,temp,n_locf);
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n_loc,n_loc,n_locf,-1.0,temp,n_locf,k_fc,n_loc,1.0,k_cc,n_loc);

	for (int i=0; i<n_loc*n_loc; i++)
		k_e[i] = k_cc[i];

	cblas_dgemv(CblasRowMajor,CblasNoTrans,n_loc,n_locf,-1.0,temp,n_loc,f_f,1,1.0,f_c,1);

	for (int i=0; i<n_loc; i++)
		f_e[i] = f_c[i];

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
