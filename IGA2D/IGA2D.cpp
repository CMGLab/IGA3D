/*
 * IGA2D.cpp
 *
 *  Created on: May 26, 2015
 *      Author: Christopher Coley
 */

#include <iostream>
#include <cmath>
using namespace std;

struct mesh {
	int n_q;
	int dim;
	int p1;
	int p2;
	int n1;
	int n2;
	int n_el_1;
	int n_el_2;
	int n_el;
	double* Xi_1;
	double* Xi_2;
	double** P;
	double* w;
	int** IEN_1;
	int** IEN_2;
	int** IEN;
	double** C_operators1;
	double** C_operators2;
	double*** C_operators;
	double*** P_b;
	double** w_b;
	double* xi_q;
	double* w_q;
	double** K;
	double* F;
	double kappa;
	double f;
	double u_R;
	double beta;
	int* BC;
	double* g;
	int** Neumann;
	double** h;
	int** Robin;

};

int findspan(int p, int n, double Xi[], double xi);
void computesplinebasis(int p, int n, double Xi[], double xi, double N[]);
int ExtractBasis1D(int p, int n, double Xi[], double** C_operators);
mesh ExtractBasis(mesh mesh1);
mesh ExtractGeometry(mesh mesh1);
void IEN1D(int p, int n, double Xi[], int** IEN);
void IEN2D(mesh mesh1);
//double NURBS2D(mesh mesh1, double xi_1, double xi_2, int i1, int i2, double N1[], double N2[]);
//void NURBScurve(mesh mesh1, double N[], double curve[][11]);
//double NURBS1D(mesh mesh1, double xi, int i, double N[]);
mesh QuadratureData(mesh mesh1);
mesh ProblemInitialization(mesh mesh1);
mesh DefineBCs(mesh mesh1);
void ElementFormation(double** k_e, double* f_e, int e, mesh mesh1);
mesh ElementAssembly(mesh mesh1, int e, double** k_e, double* f_e);
void ShapeFunction(mesh mesh1, int e, double xi_1, double xi_2, int n_loc, double* R, double** dRdx, double* x, double** J);
void BernsteinBasisAndDerivs(int p1, int p2, double xi_1, double xi_2, double* B, double** dBdxi);
void BernsteinBasisAndDeriv(int p, double xi, double* N, double* dNdxi);

int main() {

	mesh mesh1;
	double xi = 0.25;

	mesh1.dim = 2;
	mesh1.p1 = 2;
	mesh1.p2 = 2;
	mesh1.n1 = 4;
	mesh1.n2 = 4;
	mesh1.Xi_1 = new double[mesh1.n1+mesh1.p1+1];
	mesh1.Xi_2 = new double[mesh1.n2+mesh1.p2+1];
	mesh1.P = new double*[mesh1.dim];
	for (int i=0; i<mesh1.dim; i++)
		mesh1.P[i] = new double[mesh1.n1*mesh1.n2];
	mesh1.w = new double[mesh1.n1*mesh1.n2];

//	double N[mesh1.n1];
//	double curve[mesh1.dim][11] = {};

	//double Xi_init[mesh1.n1+mesh1.p1+1] = {0, 0, 0.5, 1, 1};
	double Xi_init[mesh1.n1+mesh1.p1+1] = {0, 0, 0, 0.5, 1, 1, 1};
	//double w_init[mesh1.n1*mesh1.n2] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
	double w_init[mesh1.n1*mesh1.n2] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	//double P_init[mesh1.dim][mesh1.n1*mesh1.n2] = {{0,0,0,0.5,0.5,0.5,1,1,1},{0,0.5,1,0,0.5,1,0,0.5,1}};
	double P_init[mesh1.dim][mesh1.n1*mesh1.n2] = {{0,0,0,0,0.3,0.3,0.3,0.3,0.6,0.6,0.6,0.6,1,1,1,1},{0,0.3,0.6,1,0,0.3,0.6,1,0,0.3,0.6,1,0,0.3,0.6,1}};

	for (int i=0; i<mesh1.n1+mesh1.p1+1; i++)
		mesh1.Xi_1[i] = Xi_init[i];

	for (int i=0; i<mesh1.n2+mesh1.p2+1; i++)
		mesh1.Xi_2[i] = Xi_init[i];

	for (int i=0; i<mesh1.n1*mesh1.n2; i++)
			mesh1.w[i] = w_init[i];

	for (int i=0; i<mesh1.dim; i++){
		for (int j=0; j<mesh1.n1*mesh1.n2; j++){
			mesh1.P[i][j] = P_init[i][j];
		}
	}

	//double xi_1 = 0.25;
	//double xi_2 = 0.25;
	//int i1 = 1;
	//int i2 = 1;
	//double N1[mesh1.n1];
	//double N2[mesh1.n2];
	//double R;

	//R = NURBS2D(mesh1,xi_1,xi_2,i1,i2,N1,N2);

	mesh1 = ProblemInitialization(mesh1);
	mesh1 = QuadratureData(mesh1);

	int counter = 1;
	int temp = mesh1.Xi_1[0];
	for (int i=0; i<mesh1.n1+mesh1.p1+1; i++) {
		if (temp != mesh1.Xi_1[i]) {
			counter++;
			temp = mesh1.Xi_1[i];
		}
	}

	mesh1.C_operators1 = new double*[mesh1.p1+1];
	for (int i=0; i<mesh1.p1+1; i++)
		mesh1.C_operators1[i] = new double[(mesh1.p1+1)*(counter-1)];

	counter = 1;
	temp = mesh1.Xi_2[0];
	for (int i=0; i<mesh1.n2+mesh1.p2+1; i++) {
		if (temp != mesh1.Xi_2[i]) {
			counter++;
			temp = mesh1.Xi_2[i];
		}
	}

	mesh1.C_operators2 = new double*[mesh1.p2+1];
	for (int i=0; i<mesh1.p2+1; i++)
		mesh1.C_operators2[i] = new double[(mesh1.p2+1)*(counter-1)];

	mesh1 = ExtractBasis(mesh1);

	mesh1.IEN_1 = new int*[mesh1.p1+1];
	for (int i=0; i<mesh1.p1+1; i++)
		mesh1.IEN_1[i] = new int[mesh1.n_el_1];

	mesh1.IEN_2 = new int*[mesh1.p2+1];
	for (int i=0; i<mesh1.p2+1; i++)
		mesh1.IEN_2[i] = new int[mesh1.n_el_2];

	mesh1.IEN = new int*[mesh1.p1*mesh1.p2+mesh1.p1+mesh1.p2+1];
	for (int i=0; i<mesh1.p1*mesh1.p2+mesh1.p1+mesh1.p2+1; i++)
		mesh1.IEN[i] = new int[mesh1.n_el];

	IEN2D(mesh1);

	mesh1 = ExtractGeometry(mesh1);

	mesh1 = DefineBCs(mesh1);

	mesh1.K = new double* [mesh1.n1*mesh1.n2];
	for (int i=0; i<mesh1.n1*mesh1.n2; i++)
		mesh1.K[i] = new double [mesh1.n1*mesh1.n2];
	for (int i=0; i<mesh1.n1*mesh1.n2; i++) {
		for (int j=0; j<mesh1.n1*mesh1.n2; j++)
			mesh1.K[i][j] = 0;
	}

	mesh1.F = new double [mesh1.n1*mesh1.n2];
	for (int i=0; i<mesh1.n1*mesh1.n2; i++)
		mesh1.F[i] = 0;

	double** k_e;
	k_e = new double* [(mesh1.p1+1)*(mesh1.p2+1)];
	for (int i=0; i<(mesh1.p1+1)*(mesh1.p2+1); i++)
		k_e[i] = new double [(mesh1.p1+1)*(mesh1.p2+1)];

	double* f_e;
	f_e = new double [(mesh1.p1+1)*(mesh1.p2+1)];

	for (int e=0; e<mesh1.n_el; e++) {
		ElementFormation(k_e, f_e, e, mesh1);
		mesh1 = ElementAssembly(mesh1, e, k_e, f_e);
	}

	for (int i=0; i<mesh1.n1*mesh1.n2; i++) {
		if (mesh1.BC[i] == 1) {
			mesh1.K[i][i] = 1;
			mesh1.F[i] = mesh1.g[i];
		}
	}

	/*for (int i=0; i<mesh1.n1*mesh1.n2; i++) {
		for (int j=0; j<mesh1.n1*mesh1.n2; j++) {
			cout << mesh1.K[i][j] << ", ";
		}
		cout << endl;
	}

	cout << "--------" << endl;

	for (int i=0; i<mesh1.n1*mesh1.n2; i++)
		cout << mesh1.F[i] << ", ";
	cout << endl;*/

	cout << "Stiffness matrix and forcing vector calculated" << endl;
	return 0;
}

int findspan(int p, int n, double Xi[], double xi) {
	/* Determine the knot span index */
	/* Based on Algorithm A2.1 of the NURBS book */

	int low, mid, high;

	/* Special case if xi is found at the end of Xi */
	if (xi == Xi[n+1])
		return n;

	/* Perform binary search */
	low = p;
	high = n+1;
	mid = (low + high)/2;
	while (xi < Xi[mid] || xi >= Xi[mid+1]) {
		if (xi < Xi[mid])
			high = mid;
		else
			low = mid;
		mid = (low + high)/2;
	}

	return mid;
}

void computesplinebasis(int p, int n, double Xi[], double xi, double N[]) {
	/* Compute basis functions */
	/* Based on Algorithm A2.2 of the NURBS book */

	int i,j,r;
	int l;
	double saved, temp;
	double left[p+1], right[p+1];
	double Ntemp[p];

	/* Reset N */
	for (i=0; i<=n-1; i++)
		N[i] = 0;

	/* Find knot span index such that u_l <= xi < u_(l+1)*/

	l = findspan(p,n,Xi,xi);

	/* If xi is at end of Xi, simply set basis function to 1 */

	if (l == n) {
		N[n-1] = 1;
		return;
	}

	/* Otherwise, calculate non-zero basis functions */
	Ntemp[0] = 1.0;

	for (j=1; j<=p; j++) {
		left[j] = xi - Xi[l+1-j];
		right[j] = Xi[l+j] - xi;
		saved = 0.0;
		for (r=0; r<j; r++) {
			temp = Ntemp[r]/(right[r+1]+left[j-r]);
			Ntemp[r] = saved + right[r+1]*temp;
			saved = left[j-r]*temp;
		}
		Ntemp[j] = saved;
	}

	/* Insert non-zero basis functions into array of all basis functions */
	i = l-p;
	for (j=0; j<=p; j++) {
		N[i] = Ntemp[j];
		i++;
	}
}

int ExtractBasis1D(int p, int n, double Xi[], double** C_operators) {
	/* Constructs the element extraction operators and corresponding IEN array for a one-dimensional B-spline basis */

	int a = p+1;
	int b = a+1;
	int n_el = 1;
	int i, j;
	int k, mult, r, s, save;
	double numer, a1;
	double alphas[p];
	double** C_current = new double*[p+1];
	for (int i=0; i<p+1; i++)
		C_current[i] = new double[p+1];

	for (i=0; i<p+1; i++) {
		for (j=0; j<p+1; j++) {
			if (i == j)
				C_current[i][j] = 1;
			else
				C_current[i][j] = 0;
		}
	}

	double** C_next = new double*[p+1];
	for (i=0; i<p+1; i++)
		C_next[i] = new double[p+1];

//	int counter = 1;
//	int temp = Xi[0];
//	for (i=0; i<n+p+1; i++) {
//		if (temp != Xi[i]) {
//			counter++;
//			temp = Xi[i];
//		}
//	}
//
//	C_operators = new double*[p+1];
//	for (i=0; i<p+1; i++)
//		C_operators[i] = new double[(p+1)*(counter-1)];

	while (b < n+p+1) {
		for (i=0; i<p+1; i++) {
			for (j=0; j<p+1; j++) {
				if (i == j)
					C_next[i][j] = 1;
				else
					C_next[i][j] = 0;
			}
		}

		i = b;

		while (b < n+p+1 && Xi[b] == Xi[b-1])
			b = b+1;

		mult = b-i+1;

		if (mult < p) {
			numer = Xi[b-1] - Xi[a-1];

			for (j=p; j>=mult+1; j--) {

				alphas[j-mult-1] = numer/(Xi[a+j-1] - Xi[a-1]);
			}

			r = p - mult;

			for (j=1; j<=r; j++) {
				save = r-j+1;
				s =  mult+j;

				for (int k=p+1; k>=s+1; k--) {
					a1 = alphas[k-s-1];
					for (int l=0; l<p+1; l++)
						C_current[l][k-1] = a1*C_current[l][k-1] + (1-a1)*C_current[l][k-2];
				}

				if (b < n+p+1)
					for (int l=0; l<=j; l++) {
						C_next[save+l-1][save-1] = C_current[p-j+l][p];
					}
			}
		}

		for (i=0; i<p+1; i++) {
			for (j=0; j<p+1; j++)
				C_operators[i][j+(n_el-1)*(p+1)] = C_current[i][j];
		}

		for (i=0; i<p+1; i++) {
			for (j=0; j<p+1; j++)
				C_current[i][j] = C_next[i][j];
		}

		if (b < n+p+1) {
			a = b;
			b = a+1;
			n_el++;
		}
	}

	return n_el;
}

mesh ExtractBasis(mesh mesh1) {
	/* Constructs element extraction operators for a two-dimensional B-spline basis */

	int e,e1,e2,i,j;

	mesh1.n_el_1 = ExtractBasis1D(mesh1.p1, mesh1.n1, mesh1.Xi_1, mesh1.C_operators1);
	mesh1.n_el_2 = ExtractBasis1D(mesh1.p2, mesh1.n2, mesh1.Xi_2, mesh1.C_operators2);
	mesh1.n_el = mesh1.n_el_1*mesh1.n_el_2;

	mesh1.C_operators = new double** [(mesh1.p1+1)*(mesh1.p2+1)];
	for (i=0; i<(mesh1.p1+1)*(mesh1.p2+1); i++) {
		mesh1.C_operators[i] = new double* [(mesh1.p1+1)*(mesh1.p2+1)];
		for (j=0; j<(mesh1.p1+1)*(mesh1.p2+1); j++) {
			mesh1.C_operators[i][j] = new double [mesh1.n_el];
		}
	}

	for (e1=0; e1<mesh1.n_el_1; e1++) {
		for (e2=0; e2<mesh1.n_el_2; e2++) {
			e = e1*mesh1.n_el_2 + e2;

			for (int i1=0; i1<mesh1.p1+1; i1++) {
				for (int j1=0; j1<mesh1.p1+1; j1++) {
					for (int i2=0; i2<mesh1.p2+1; i2++) {
						for (int j2=0; j2<mesh1.p2+1; j2++) {
								i = i1*(mesh1.p2+1) + i2;
								j = j1*(mesh1.p2+1) + j2;
								mesh1.C_operators[i][j][e] = mesh1.C_operators1[i1][e1*(mesh1.p1+1)+j1]*mesh1.C_operators2[i2][e2*(mesh1.p2+1)+j2];
						}
					}
				}
			}
		}
	}

	return mesh1;
}

mesh ExtractGeometry(mesh mesh1) {
	/* Computes the Bezier control points and weights corresponding to a NURBS surface */

	int i;
	double** Pw;
	Pw = new double* [mesh1.dim+1];
	for (i=0; i<mesh1.dim+1; i++)
		Pw[i] = new double[mesh1.n1*mesh1.n2];
	double P_e[(mesh1.p1+1)*(mesh1.p2+1)];
	double** C_e;
	C_e = new double* [(mesh1.p1+1)*(mesh1.p2+1)];
	for (int i=0; i<(mesh1.p1+1)*(mesh1.p2+1); i++)
		C_e[i] = new double[(mesh1.p1+1)*(mesh1.p2+1)];
	double*** P_bw;
	P_bw = new double** [mesh1.dim+1];
	for (int i=0; i<mesh1.dim+1; i++) {
		P_bw[i] = new double* [(mesh1.p1+1)*(mesh1.p2+1)];
		for (int j=0; j<(mesh1.p1+1)*(mesh1.p2+1); j++)
			P_bw[i][j] = new double [mesh1.n_el];
	}
	mesh1.w_b = new double* [(mesh1.p1+1)*(mesh1.p2+1)];
	for (int i=0; i<(mesh1.p1+1)*(mesh1.p2+1); i++)
		mesh1.w_b[i] = new double [mesh1.n_el];
	mesh1.P_b = new double** [mesh1.dim];
	for (int i=0; i<mesh1.dim; i++) {
		mesh1.P_b[i] = new double* [(mesh1.p1+1)*(mesh1.p2+1)];
		for (int j=0; j<(mesh1.p1+1)*(mesh1.p2+1); j++)
			mesh1.P_b[i][j] = new double [mesh1.n_el];
	}

	for (i=0; i<mesh1.dim; i++) {
		for (int j=0; j<mesh1.n1*mesh1.n2; j++) {
			Pw[i][j] = mesh1.P[i][j]*mesh1.w[j];
		}
	}

	for (int j=0; j<mesh1.n1*mesh1.n2; j++)
		Pw[mesh1.dim][j] = mesh1.w[j];

	for (int e=0; e<mesh1.n_el; e++) {
		for (int j=0; j<mesh1.dim+1; j++) {
			for (int a=0; a<(mesh1.p1+1)*(mesh1.p2+1); a++) {
				i = mesh1.IEN[a][e];
				P_e[a] = Pw[j][i-1];
			}

			for (int k=0; k<(mesh1.p1+1)*(mesh1.p2+1); k++) {
				for (int l=0; l<(mesh1.p1+1)*(mesh1.p2+1); l++)
					C_e[k][l] = mesh1.C_operators[k][l][e];
			}

			for (int k=0; k<(mesh1.p1+1)*(mesh1.p2+1); k++) {
				P_bw[j][k][e] = 0;
				for (int l=0; l<(mesh1.p1+1)*(mesh1.p2+1); l++)
					P_bw[j][k][e] += C_e[l][k] * P_e[l];
			}
		}
	}

	//divide out weights and extract new control points and weights
	for (int e=0; e<mesh1.n_el; e++) {
		for (int i=0; i<(mesh1.p1+1)*(mesh1.p2+1); i++)
			mesh1.w_b[i][e] = P_bw[mesh1.dim][i][e];
		for (int j=0; j<mesh1.dim; j++) {
			for (int k=0; k<(mesh1.p1+1)*(mesh1.p2+1); k++)
				mesh1.P_b[j][k][e] = P_bw[j][k][e]/mesh1.w_b[k][e];
		}
	}

	return mesh1;
}

void IEN1D(int p, int n, double Xi[], int** IEN) {
	/* Computes the IEN array for a one-dimensional B-spline basis */

	int l = p+1;
	int e = 1;

	while (l<n+1) {
		for (int a=1; a<=p+1; a++) {
			IEN[a-1][e-1] = (l+a) - (p+1);
		}

		l = l+1;

		while (Xi[l] == Xi[l-1] && l < n+1)
			l++;

		if (l<n+1)
			e++;
	}
}

void IEN2D(mesh mesh1) {
	/* Computes the IEN array for a two-dimensional B-spline basis */

	int i, i1, i2, e, a;

	IEN1D(mesh1.p1, mesh1.n1, mesh1.Xi_1, mesh1.IEN_1);
	IEN1D(mesh1.p2, mesh1.n2, mesh1.Xi_2, mesh1.IEN_2);

	for (int e1=0; e1<mesh1.n_el_1; e1++){
		for (int a1=0; a1<mesh1.p1+1; a1++){

			i1 = mesh1.IEN_1[a1][e1];

			for (int e2=0; e2<mesh1.n_el_2; e2++){
				for (int a2=0; a2<mesh1.p2+1; a2++){

					i2 = mesh1.IEN_2[a2][e2];
					e = e1*mesh1.n_el_2 + e2;
					a = a1*(mesh1.p2+1) + a2;
					i = (i1-1)*mesh1.n2 + i2;

					mesh1.IEN[a][e] = i;
				}
			}
		}
	}
}

mesh QuadratureData(mesh mesh1) {

	if (mesh1.p1 > mesh1.p2)
		mesh1.n_q = mesh1.p1+1;
	else
		mesh1.n_q = mesh1.p2+1;

	switch (mesh1.n_q) {
	case 1:
		mesh1.xi_q = new double [1];
		mesh1.w_q = new double [1];
		mesh1.xi_q[0] = 1./2.;
		mesh1.w_q[0] = 1.;
		break;
	case 2:
		mesh1.xi_q = new double [2];
		mesh1.w_q = new double [2];
		mesh1.xi_q[0] = 1./2.-1./2.*sqrt(1./3.);
		mesh1.xi_q[1] = 1./2.+1./2.*sqrt(1./3.);
		mesh1.w_q[0] = 1./2.;
		mesh1.w_q[1] = 1./2.;
		break;
	case 3:
		mesh1.xi_q = new double [3];
		mesh1.w_q = new double [3];
		mesh1.xi_q[0] = 1./2.-1./2.*sqrt(3./5.);
		mesh1.xi_q[1] = 1./2.;
		mesh1.xi_q[2] = 1./2.+1./2.*sqrt(3./5.);
		mesh1.w_q[0] = 5./18.;
		mesh1.w_q[1] = 4./9.;
		mesh1.w_q[2] = 5./18.;
		break;
	case 4:
		mesh1.xi_q = new double [4];
		mesh1.w_q = new double [4];
		mesh1.xi_q[0] = 1./2.-1./2.*sqrt(3./7.-2./7.*sqrt(6./5.));
		mesh1.xi_q[1] = 1./.2+1./2.*sqrt(3./7.-2./7.*sqrt(6./5.));
		mesh1.xi_q[2] = 1./2.-1./2.*sqrt(3./7.+2./7.*sqrt(6./5.));
		mesh1.xi_q[3] = 1./2.+1./2.*sqrt(3./7.+2./7.*sqrt(6./5.));
		mesh1.w_q[0] = 1./2.*(18.+sqrt(30.))/36.;
		mesh1.w_q[1] = 1./2.*(18.+sqrt(30.))/36.;
		mesh1.w_q[2] = 1./2.*(18.-sqrt(30.))/36.;
		mesh1.w_q[3] = 1./2.*(18.-sqrt(30.))/36.;
		break;
	}

	return mesh1;
}

/*void NURBScurve(mesh mesh1, double N[], double curve[][11]) {
	// Generates a NURBS curve

	int d = mesh1.dim;
	int p = mesh1.p;
	int n = mesh1.n;
	double Xi[n+p+1];
	for (int i=0; i<n+p+1; i++)
		Xi[i] = mesh1.Xi[i];
	double P[d][n];
	for (int i=0; i<d; i++) {
		for (int j=0; j<n; j++) {
			P[i][j] = mesh1.P[i][j];
		}
	}
	double w[n];
	for (int i=0; i<n; i++)
		w[i] = mesh1.w[i];

	int i,j,k;
	int lenxi;
	double dt = 0.1;
	double R[n] = {};
	double sum;
	double Xistart, Xiend;

	// Reset curve
	for (j=0; j<=2; j++) {
		for (i=0; i<=d-1; i++) {
			curve[i][j] = 0;
		}
	}

	Xistart = Xi[0];
	Xiend = Xi[n+p];

	lenxi = (Xiend - Xistart)/dt;

	double xi[lenxi];

	for (j=0; j<=lenxi; j++) {
		xi[j] = Xistart + (Xiend - Xistart)*j/lenxi;
	}

	for (j=0; j<=lenxi; j++) {
		for (i=0; i<=n-1; i++)
			R[i] = NURBS1D(mesh1,xi[j],i,N);
		for (i=0; i<=d-1; i++) {
			sum = 0.0;
			for (k=0; k<=n-1; k++)
				sum += P[i][k]*R[k];
			curve[i][j] = sum;
		}
	}
}*/

/*double NURBS1D(mesh mesh1, double xi, int i, double N[]) {
	// Evaluates one-dimensional NURBS basis functions at specified parameter points

	int p = mesh1.p;
	int n = mesh1.n;
	double Xi[n+p+1];
	for (int k=0; k<n+p+1; k++)
		Xi[k] = mesh1.Xi[k];
	double w[n];
	for (int k=0; k<n; k++)
		w[k] = mesh1.w[k];

	int j;
	double R;
	double sum = 0.0;

	computesplinebasis(mesh1,xi,N);

	for (j=0; j<=n; j++) {
		sum += w[j]*N[j];
	}

	R = w[i]*N[i]/sum;

	return R;
}*/

/*double NURBS2D(mesh mesh1, double xi_1, double xi_2, int i1, int i2, double N1[], double N2[]) {
	// Evaluates two-dimensional NURBS basis functions at specified parameter points

	int p1 = mesh1.p1;
	int p2 = mesh1.p2;
	int n1 = mesh1.n1;
	int n2 = mesh1.n2;
	double Xi_1[n1+p1+1];
	for (int k=0; k<n1+p1+1; k++)
		Xi_1[k] = mesh1.Xi_1[k];
	double Xi_2[n2+p2+1];
		for (int k=0; k<n2+p2+1; k++)
			Xi_2[k] = mesh1.Xi_2[k];
	double w[n1*n2];
	for (int k=0; k<n1*n2; k++)
		w[k] = mesh1.w[k];

	int i,j;
	double R;
	double sum = 0.0;
	double N_temp[n1*n2];

	computesplinebasis(p1,n1,Xi_1,xi_1,N1);
	computesplinebasis(p2,n2,Xi_2,xi_2,N2);

	//kronecker product
	for (j=0; j<n2; j++){
		for (i=0; i<n1; i++){
			N_temp[i+n1*j] = N1[i]*N2[j];
		}
	}

	for (j=0; j<=n1*n2; j++) {
		sum += w[j]*N_temp[j];
	}

	R = w[i1+n1*i2]*N_temp[i1+n1*i2]/sum;

	return R;
}*/

mesh ProblemInitialization(mesh mesh1) {
	mesh1.kappa = 385;
	mesh1.f = 0;
	mesh1.u_R = 0;
	mesh1.beta = 0;

	return mesh1;
}

mesh DefineBCs(mesh mesh1) {
	int j,k;
	double u_o = 70.;
	double u_i = 200.;

	mesh1.BC = new int [mesh1.n1*mesh1.n2];
	for (int i=0; i<mesh1.n1*mesh1.n2; i++)
		mesh1.BC[i] = 0;

	mesh1.g = new double [mesh1.n1*mesh1.n2];
	for (int i=0; i<mesh1.n1*mesh1.n2; i++)
			mesh1.g[i] = 0;

	for (int i=0; i<mesh1.n1; i++) {
		j = i*mesh1.n2;
		k = i*mesh1.n2 + mesh1.n2 - 1;

		mesh1.g[j] = u_i;
		mesh1.g[k] = u_o;

		mesh1.BC[j] = 1;
		mesh1.BC[k] = 1;
	}

	mesh1.Neumann = new int* [4];
	for (int i=0; i<4; i++) {
		mesh1.Neumann[i] = new int [mesh1.n_el];
	}
	for (int i=0; i<4; i++) {
		for (int j=0; j<mesh1.n_el; j++)
			mesh1.Neumann[i][j] = 0;
	}

	mesh1.h = new double* [4];
		for (int i=0; i<4; i++)
			mesh1.h[i] = new double [mesh1.n_el];
	for (int i=0; i<4; i++) {
		for (int j=0; j<mesh1.n_el; j++)
			mesh1.h[i][j] = 0;
	}

	for (int i=0; i<mesh1.n_el_2; i++)
		mesh1.Neumann[3][i] = 1;

	for (int i=mesh1.n_el-1; i>mesh1.n_el-mesh1.n_el_2; i--)
		mesh1.Neumann[1][i] = 1;

	mesh1.Robin = new int* [4];
		for (int i=0; i<4; i++)
			mesh1.Robin[i] = new int [mesh1.n_el];
	for (int i=0; i<4; i++) {
		for (int j=0; j<mesh1.n_el; j++)
			mesh1.Robin[i][j] = 0;
	}

	return mesh1;
}

void ElementFormation(double** k_e, double* f_e, int e, mesh mesh1) {
	int n_loc = (mesh1.p1+1)*(mesh1.p2+1);

	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<n_loc; j++)
		k_e[i][j] = 0;
	}

	for (int i=0; i<n_loc; i++)
	f_e[i] = 0;

	double* R_e;
	R_e = new double [n_loc];
	double** dRdx;
	dRdx = new double* [n_loc];
	for (int i=0; i<n_loc; i++)
		dRdx[i] = new double [2];

	double** J;
	J = new double* [2];
	for (int i=0; i<2; i++)
		J[i] = new double [2];
	double* x;
	x = new double [2];
	x[0] = 0;
	x[1] = 0;

	double xi_t1;
	double xi_t2;
	double t_t[2];

	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<n_loc; j++)
			k_e[i][j] = 0;
	}
	for (int i=0; i<n_loc; i++)
		f_e[i] = 0;

	for (int q1=0; q1<mesh1.n_q; q1++) {
		for (int q2=0; q2<mesh1.n_q; q2++) {

			ShapeFunction(mesh1, e, mesh1.xi_q[q1], mesh1.xi_q[q2], n_loc, R_e, dRdx, x, J);

			double j = J[0][0]*J[1][1] - J[0][1]*J[1][0];
			double k_loc = mesh1.kappa;
			double f_loc = mesh1.f;

			for (int a=0; a<n_loc; a++) {
				for (int b=0; b<n_loc; b++) {
					k_e[a][b] += k_loc*(dRdx[a][0]*dRdx[b][0] + dRdx[a][1]*dRdx[b][1])*mesh1.w_q[q1]*mesh1.w_q[q2]*j;
				}
				f_e[a] += R_e[a]*f_loc*mesh1.w_q[q1]*mesh1.w_q[q2]*j;
			}
		}
	}

	for (int side=0; side<4; side++) {
		if (mesh1.Neumann[side][e] == 1 || mesh1.Robin[side][e] == 1) {
			for (int q=0; q<mesh1.n_q; q++) {
				switch (side) {
				case 1:
					xi_t1 = mesh1.xi_q[q];
					xi_t2 = 0;
					t_t[0] = 1;
					t_t[1] = 0;
					break;
				case 2:
					xi_t1 = 1;
					xi_t2 = mesh1.xi_q[q];
					t_t[0] = 0;
					t_t[1] = 1;
					break;
				case 3:
					xi_t1 = mesh1.xi_q[q];
					xi_t2 = 1;
					t_t[0] = 1;
					t_t[1] = 0;
					break;
				case 4:
					xi_t1 = 0;
					xi_t2 = mesh1.xi_q[q];
					t_t[0] = 0;
					t_t[1] = 1;
					break;
				}

				ShapeFunction(mesh1, e, xi_t1, xi_t2, n_loc, R_e, dRdx, x, J);
				double j_prod[2] = {0,0};
				for (int i=0; i<2; i++) {
					for (int j=0; j<2; j++)
						j_prod[i] += J[i][j]*t_t[j];
				}
				double j_boundary = sqrt(j_prod[0]*j_prod[0]+j_prod[1]*j_prod[1]);

				if (mesh1.Neumann[side][e] == 1) {
					double h_loc = mesh1.h[side][e];
					for (int a=0; a<n_loc; a++)
						f_e[a] += R_e[a]*h_loc*mesh1.w_q[q]*j_boundary;
				} else if (mesh1.Robin[side][e] == 1) {
					double beta_loc = mesh1.beta;
					double u_r_loc = mesh1.u_R;
					for (int a=0; a<n_loc; a++) {
						for (int b=0; b<n_loc; b++)
							k_e[a][b] += beta_loc*R_e[a]*R_e[b]*mesh1.w_q[q]*j_boundary;
					}
				}
			}
		}
	}
}

mesh ElementAssembly(mesh mesh1, int e, double** k_e, double* f_e) {
	int n_loc = (mesh1.p1+1)*(mesh1.p2+1);

	for (int a=0; a<n_loc; a++) {
		int i = mesh1.IEN[a][e] - 1;
		if (mesh1.BC[i] == 0) {
			for (int b=0; b<n_loc; b++) {
				int j = mesh1.IEN[b][e] - 1;
				if (mesh1.BC[j] == 0){
					mesh1.K[i][j] += k_e[a][b];
				} else {
					mesh1.F[i] -= k_e[a][b]*mesh1.g[j];
				}
			}
			mesh1.F[i] += f_e[a];
		}
	}

	return mesh1;
}

void ShapeFunction(mesh mesh1, int e, double xi_1, double xi_2, int n_loc, double* R, double** dRdx, double* x, double** J) {

	for (int i=0; i<n_loc; i++)
		R[i] = 0;
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<2; j++)
			dRdx[i][j] = 0;
	}

	double wb = 0;
	double dwbdxi[2] = {0, 0};
	double** dRdxi;
	dRdxi = new double* [n_loc];
	for (int i=0; i<n_loc; i++)
		dRdxi[i] = new double [2];
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<2; j++)
			dRdxi[i][j] = 0;
	}

	double** dxdxi;
	dxdxi = new double* [2];
	for (int i=0; i<2; i++)
		dxdxi[i] = new double [2];
	for (int i=0; i<2; i++) {
		for (int j=0; j<2; j++)
			dxdxi[i][j] = 0;
	}

	double** dxidx;
	dxidx = new double* [2];
	for (int i=0; i<2; i++)
		dxidx[i] = new double [2];
	for (int i=0; i<2; i++) {
		for (int j=0; j<2; j++)
			dxidx[i][j] = 0;
	}

	double* B;
	B = new double [n_loc];
	double** dBdxi;
	dBdxi = new double* [n_loc];
	for (int i=0; i<n_loc; i++)
		dBdxi[i] = new double [2];

	BernsteinBasisAndDerivs(mesh1.p1, mesh1.p2, xi_1, xi_2, B, dBdxi);

	//weighting and derivatives
	for (int a=0; a<n_loc; a++){
		wb += mesh1.w_b[a][e]*B[a];
		dwbdxi[0] += mesh1.w_b[a][e]*dBdxi[a][0];
		dwbdxi[1] += mesh1.w_b[a][e]*dBdxi[a][1];
	}

	//basis functions and parametric derivatives
	for (int a=0; a<n_loc; a++) {
		for (int b=0; b<n_loc; b++) {
			R[a] += mesh1.w[mesh1.IEN[a][e]-1] * mesh1.C_operators[a][b][e]*B[b]/wb;
			for (int i=0; i<2; i++)
				dRdxi[a][i] += mesh1.w[mesh1.IEN[a][e]-1] * mesh1.C_operators[a][b][e] * (dBdxi[b][i]/wb - dwbdxi[i]*B[b]/(wb*wb));
		}
	}

	//physical space quantities
	for (int a=0; a<n_loc; a++) {
		for (int i=0; i<2; i++)
			x[i] += mesh1.w_b[a][e]*mesh1.P_b[i][a][e]*B[a]/wb;
		for (int i=0; i<2; i++) {
			for (int j=0; j<2; j++)
				dxdxi[i][j] += mesh1.w_b[a][e]*mesh1.P_b[i][a][e]*(dBdxi[a][j]/wb - dwbdxi[j]*B[a]/(wb*wb));
		}
	}

	for (int i=0; i<2; i++) {
		for (int j=0; j<2; j++)
			J[i][j] = dxdxi[i][j];
	}

	double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
	dxidx[0][0] = J[1][1]/detJ;
	dxidx[0][1] = -J[0][1]/detJ;
	dxidx[1][0] = -J[1][0]/detJ;
	dxidx[1][1] = J[0][0]/detJ;

	//physical derivatives
	for (int a=0; a<n_loc; a++) {
		for (int i=0; i<2; i++) {
			for (int j=0; j<2; j++)
				dRdx[a][i] += dRdxi[a][j]*dxidx[j][i];
		}
	}
}

void BernsteinBasisAndDerivs(int p1, int p2, double xi_1, double xi_2, double* B, double** dBdxi) {
	double* N;
	N = new double [p1+1];
	double* dNdxi;
	dNdxi = new double [p1+1];
	double* M;
	M = new double [p2+1];
	double* dMdxi;
	dMdxi = new double [p2+1];

	BernsteinBasisAndDeriv(p1, xi_1, N, dNdxi);
	BernsteinBasisAndDeriv(p2, xi_2, M, dMdxi);

	int a = 0;

	for (int i=0; i<p1+1; i++) {
		for (int j=0; j<p2+1; j++) {
			B[a] = N[i]*M[j];
			dBdxi[a][0] = dNdxi[i]*M[j];
			dBdxi[a][1] = N[i]*dMdxi[j];
			a++;
		}
	}
}

void BernsteinBasisAndDeriv(int p, double xi, double* N, double* dNdxi) {
	int n = 1; //can change this later to get higher derivatives

	double ndu[p+1][p+1];
	ndu[0][0] = 1;
	double ders[2][p+1];

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
			r *= (p-k);
		}
	}

	for (int j=0; j<=p; j++) {
		N[j] = ders[0][j];
		dNdxi[j] = ders[1][j];
	}
}
