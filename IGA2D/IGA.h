/*
 * IGA.h
 *
 *  Created on: Aug 19, 2015
 *      Author: Christopher
 */

#ifndef IGA_H_
#define IGA_H_

struct ASG_Element {
	int p;
	//int pf;
	int nloc;
	int elemType;
	double* points;
	double* weights;
	double* bezierPoints;
	double* bezierWeights;
	double* extraction;
	int* connectivity;
};

struct BC_Struct {
	int* Neumann;
	double* h;
	int* Robin;
	double* beta;
};

struct Prob_Params {
	double* kappa;
	double* kappaArt;
	double* f;
	double* u;
	double C;
};

typedef struct _p_Vec* Vec;
typedef struct _p_Mat* Mat;

void Extract_Geometry(ASG_Element* mesh, int dim, int numElems, double** globalBezierPoints, double** globalBezierWeights);
void Quadrature_Data(int numQ, double* quadPts, double* quadWts);
void Bernstein_Basis_and_Deriv(double** ders, int p, double xi, int n);
void Bernstein_Basis_and_Derivs(double** B, int p1, int p2, double xi_1, double xi_2, int n);
void Full_Bernstein_Basis_and_Derivs(double*** fullBasis, int nloc, int p, int n_q, int n, double* quadPts);
void Boundary_Bernstein_Basis_and_Derivs(double*** fullBasis, int nloc, int p, int n_q, int n, double* quadPts);
void ShapeFunction(ASG_Element element, double basis[], double R[], double dRdx[], double d2Rdx2[], double x[], double J[], int flag);
double Grid_Spacing(double bezierPoints[], int p);
void Element_Assembly(ASG_Element element, int* BC, double* g, double* k, double* f, Mat K, Vec F);
double Get_Local_Value(double vi[], double basis[], int nloc);
void Extract_Concentration(ASG_Element* mesh, double*** fullBasis, int numElems, double* plotPts, int numPlotPts, double* d, double* X, double* Y, double* U);
void Subelement_Connectivity(int* subConnectivity, int numPlotPts, int numSubElems);

#endif /* IGA_H_ */
