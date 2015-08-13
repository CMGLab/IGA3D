/*
 * AD2D.cpp
 *
 *  Created on: July 3, 2015
 *  Last modified on: 12 August, 2015
 *      Author: Christopher Coley
 */

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <petscksp.h>
#include <cblas.h>
#include <lapacke.h>
using namespace std;

struct ASG_Element {
	int* p;
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

void Read_IGA_File(ASG_Element*& mesh, int& dim, int& numElems, int& numNodes);
void Read_BC_File(BC_Struct* BCs, int* BC, double* g, int numNodes, int numElems);
void Read_Problem_Parameters(Prob_Params* probParams, ASG_Element* mesh, int numNodes, int numElems, double pf);
void Extract_Geometry(ASG_Element* mesh, int dim, int numElems);
void Quadrature_Data(int numQ, double* quadPts, double* quadWts);
void Bernstein_Basis_and_Deriv(double** ders, int p, double xi, int n);
void Bernstein_Basis_and_Derivs(double** B, int p1, int p2, double xi_1, double xi_2, int n);
void Full_Bernstein_Basis_and_Derivs(double*** fullBasis, int nloc, int p, int n_q, int n, double* quadPts);
void Boundary_Bernstein_Basis_and_Derivs(double*** fullBasis, int nloc, int p, int n_q, int n, double* quadPts);
void ShapeFunction(ASG_Element element, double** basis, double* R, double** dRdx, double** d2Rdx2, double* x, double* J, int flag);
void ShapeFine(ASG_Element element, int pf, double** basis, double** dRfdx, double* J);
double Grid_Spacing(double bezierPoints[], int p);
void Element_Formation(ASG_Element element, BC_Struct elementBC, Prob_Params probParams, double*** fullBasis, double*** fullBasisFine, double*** fullBoundaryBasis, double*** fullBoundaryBasisFine, int n_q, double* quadWts, int pf, double* k_e, double* f_e);
void Element_Assembly(ASG_Element element, int* BC, double* g, double* k, double* f, Mat K, Vec F);
double Get_Local_Value(double vi[], double basis[], int nloc);
void Solve_Fine_Scales(ASG_Element element, BC_Struct elementBC, Prob_Params probParams, double*** fullBasis, double*** fullBasisFine, double*** fullBoundaryBasis, double*** fullBoundaryBasisFine, int n_q, double* quadWts, int pf, double* d, double* d_f);
void Extract_Concentration(ASG_Element* mesh, double*** fullBasis, int numElems, double* plotPts, int numPlotPts, double* d, double* X, double* Y, double* U);
void Extract_Fine_Scales(double*** fullBasisFine, int numElems, int numPlotPts, int pf, double* d, double* U);
void Subelement_Connectivity(int* subConnectivity, int p);
void Write_VTK_File(double* X, double* Y, double* U, double* U_f, int numElems, int numPlotPts, int* subConnectivity, int numSubElems);

int main() {

	//TODO figure out how to do this properly
	PetscInitialize(0,(char***)0,(char*)0,(char*)0);

	int dim, numElems, numNodes, numQ;
	int numDerivs = 2;
	ASG_Element* mesh;

	//TODO Include fine scale data stuff somehow
	int pf = 1;

	// Read in geometry from .iga file
	Read_IGA_File(mesh,dim,numElems,numNodes);

	// Initialize boundary condition variables
	BC_Struct boundaryConditions[numElems];
	int BC[numNodes];
	double g[numNodes];
	for (int i=0; i<numElems; i++) {
		boundaryConditions[i].Neumann = new int[4];
		boundaryConditions[i].h = new double[4];
		boundaryConditions[i].Robin = new int[4];
		boundaryConditions[i].beta = new double[4];
	}

	// Read in boundary condition data from .bc file
	Read_BC_File(boundaryConditions,BC,g,numNodes,numElems);

	// Initialize problem parameter variables
	Prob_Params probParams[numElems];
	for (int i=0; i<numElems; i++) {
		probParams[i].kappa = new double[mesh[i].nloc];
		probParams[i].kappaArt = new double[mesh[i].nloc];
		probParams[i].f = new double[mesh[i].nloc];
		probParams[i].u = new double[2*mesh[i].nloc];
	}

	// Read in problem parameter information
	Read_Problem_Parameters(probParams,mesh,numNodes,numElems,pf);

	/* test code to look at problem parameter assignments
	for (int j=0; j<numElems; j++) {
		cout << "Element " << j+1 << endl;
		cout << "Kappa" << endl;
		for (int i=0; i<mesh[j].nloc; i++)
			cout << probParams[j].kappa[i] << " ";
		cout << endl << "KappaArt" << endl;
		for (int i=0; i<mesh[j].nloc; i++)
			cout << probParams[j].kappaArt[i] << " ";
		cout << endl << "f" << endl;
		for (int i=0; i<mesh[j].nloc; i++)
			cout << probParams[j].f[i] << " ";
		cout << endl << "u_1" << endl;
		for (int i=0; i<mesh[j].nloc; i++)
			cout << probParams[j].u[i] << " ";
		cout << endl << "u_2" << endl;
		for (int i=0; i<mesh[j].nloc; i++)
			cout << probParams[j].u[mesh[j].nloc+i] << " ";
		cout << endl << "C" << endl;
		cout << probParams[j].C << " ";
	} */

	// Extract Bezier control points and weights
	Extract_Geometry(mesh,dim,numElems);

	/* test code to look at geometry assignments
	for (int i=0; i<numElems; i++) {
		cout << "Element " << i+1 << endl;
		cout << "Nloc = " << mesh[i].nloc << endl;
		cout << "p = " << *mesh[i].p << endl;
		cout << "Connectivity:";
		for (int j=0; j<mesh[0].nloc; j++) {
			cout << ' ' << mesh[i].connectivity[j];
		}
		cout << endl;
		cout << "Extraction:" << endl;
		for (int k=0; k<mesh[i].nloc; k++) {
			for (int l=0; l<mesh[i].nloc; l++)
				cout << mesh[i].extraction[k*mesh[i].nloc+l] << ", ";
			cout << endl;
		}
		cout << "Control points:" << endl;
		for (int k=0; k<mesh[i].nloc; k++) {
			for (int l=0; l<2; l++)
				cout << mesh[i].points[k*2+l] << ", ";
			cout << endl;
		}
		cout << "Control weights:" << endl;
		for (int k=0; k<mesh[i].nloc; k++)
			cout << mesh[i].weights[k] << ", ";
		cout << endl << "---------------------" << endl;
	} */

	/* test code to look at Bezier control points and weights assignments
	for (int i=0; i<numElems; i++) {
		cout << "Element " << i+1 << endl;
		cout << "Bezier control points:" << endl;
		for (int k=0; k<mesh[i].nloc; k++) {
			for (int l=0; l<2; l++)
				cout << mesh[i].bezierPoints[k*2+l] << ", ";
			cout << endl;
		}
		cout << "Bezier control weights:" << endl;
		for (int k=0; k<mesh[i].nloc; k++)
			cout << mesh[i].bezierWeights[k] << ", ";
		cout << "---------------------" << endl;
	} */

	/* test code to look at calculation of Bernstein basis and derivatives
	int numFuncs = (numDerivs+1)*(numDerivs+2)/2;
	for (int k=0; k<(mesh[0].p[0]+1)*(mesh[0].p[0]+1); k++) {
		cout << "Quadrature point " << k+1 << endl;
		for (int i=0; i<mesh[0].nloc; i++) {
			for (int j=0; j<numFuncs; j++) {
				cout << fullBasis[i][j][k] << ", ";
			}
			cout << endl;
		}
		cout << "---------------------" << endl;
	} */

	// Use n quadrature points such that n=p+1
	if (mesh[0].p[0] > pf)
		numQ = mesh[0].p[0]+1;
	else
		numQ = pf+1;

	// Create variables to hold quadrature points and weights
	double quadPts[numQ];
	double quadWts[numQ];

	// Retrieve quadrature points and weights
	Quadrature_Data(numQ,quadPts,quadWts);

	// Pre-compute Bernstein basis and derivatives; assumes all elements are of the same type and polynomial degree
	double*** fullBasis;
	fullBasis = new double** [mesh[0].nloc];
	for (int i=0; i<mesh[0].nloc; i++) {
		fullBasis[i] = new double* [6];
		for (int j=0; j<6; j++)
			fullBasis[i][j] = new double [numQ*numQ];
	}
	Full_Bernstein_Basis_and_Derivs(fullBasis,mesh[0].nloc,mesh[0].p[0],numQ,numDerivs,quadPts);

	double*** fullBasisFine;
	fullBasisFine = new double** [(pf+1)*(pf+1)];
	for (int i=0; i<(pf+1)*(pf+1); i++) {
		fullBasisFine[i] = new double* [3];
		for (int j=0; j<3; j++)
			fullBasisFine[i][j] = new double [numQ*numQ];
	}
	Full_Bernstein_Basis_and_Derivs(fullBasisFine,(pf+1)*(pf+1),pf,numQ,1,quadPts);

	// Pre-compute Bernstein basis and derivatives on boundaries; assumes all elements are of the same type and polynomial degree
	double*** fullBoundaryBasis;
	fullBoundaryBasis = new double** [mesh[0].nloc];
	for (int i=0; i<mesh[0].nloc; i++) {
		fullBoundaryBasis[i] = new double* [6];
		for (int j=0; j<6; j++)
			fullBoundaryBasis[i][j] = new double [numQ*4];
	}
	Boundary_Bernstein_Basis_and_Derivs(fullBoundaryBasis,mesh[0].nloc,mesh[0].p[0],numQ,1,quadPts);

	double*** fullBoundaryBasisFine;
	fullBoundaryBasisFine = new double** [(pf+1)*(pf+1)];
	for (int i=0; i<(pf+1)*(pf+1); i++) {
		fullBoundaryBasisFine[i] = new double* [3];
		for (int j=0; j<3; j++)
			fullBoundaryBasisFine[i][j] = new double [numQ*4];
	}
	Boundary_Bernstein_Basis_and_Derivs(fullBoundaryBasisFine,(pf+1)*(pf+1),pf,numQ,1,quadPts);

	// Initialize variables
	Mat K;
	Vec F;
	Vec d;

	MatCreate(PETSC_COMM_WORLD,&K);
	MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,numNodes,numNodes);
	MatSetFromOptions(K);
	MatSetUp(K);

	VecCreate(PETSC_COMM_WORLD,&F);
	VecSetSizes(F,PETSC_DECIDE,numNodes);
	VecSetFromOptions(F);
	VecSetUp(F);

	VecCreate(PETSC_COMM_WORLD,&d);
	VecSetSizes(d,PETSC_DECIDE,numNodes);
	VecSetFromOptions(d);
	VecSetUp(d);

	// Element formation and assembly

	for (int i=0; i<numElems; i++) {

		// Initialize variables
		double k[mesh[i].nloc*mesh[i].nloc];
		double f[mesh[i].nloc];

		// Element formation
		Element_Formation(mesh[i],boundaryConditions[i],probParams[i],fullBasis,fullBasisFine,fullBoundaryBasis,fullBoundaryBasisFine,numQ,quadWts,pf,k,f);

		/* Output for looking at elemnt stiffness matrix and forcing vector
		for (int j=0; j<mesh[i].nloc; j++){
			for (int l=0; l<mesh[i].nloc; l++) {
				cout << k[j*mesh[i].nloc+l] << ", ";
			}
			cout << endl;
		}
		cout << "-------------" << endl;

		for (int j=0; j<mesh[i].nloc; j++)
			cout << f[j] << ", ";
		cout << endl << "-------------" << endl;
		*/

		// Element assembly
		Element_Assembly(mesh[i],BC,g,k,f,K,F);
	}

	// Apply Dirichlet BCs
	for (int i=0; i<numNodes; i++) {
		if (BC[i] == 1) {
			MatSetValue(K,i,i,1.0,INSERT_VALUES);
			VecSetValue(F,i,g[i],INSERT_VALUES);
		}
	}

	MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);

	VecAssemblyBegin(F);
	VecAssemblyEnd(F);

	/* Code to output system stiffness matrix and forcing vector
	MatView(K,PETSC_VIEWER_STDOUT_SELF);
	VecView(F,PETSC_VIEWER_STDOUT_SELF);
	*/

	// Solve System
	KSP ksp;

	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetOperators(ksp,K,K);
	KSPSetFromOptions(ksp);
	KSPSolve(ksp,F,d);

	/* Code to output control variables
	VecView(d,PETSC_VIEWER_STDOUT_SELF);
	*/

	// Extract values of control variables into an array of doubles
	double ctrVars[numNodes];
	int indices[numNodes];
	for (int i=0; i<numNodes; i++)
		indices[i] = i;

	VecGetValues(d,numNodes,indices,ctrVars);

	// Create variables to hold fine scale solution
	double d_f[(pf+1)*(pf+1)];
	double ctrVarsFine[numElems*(pf+1)*(pf+1)];

	// Solve for the fine scales
	for (int i=0; i<numElems; i++) {
		Solve_Fine_Scales(mesh[i],boundaryConditions[i],probParams[i],fullBasis,fullBasisFine,fullBoundaryBasis,fullBoundaryBasisFine,numQ,quadWts,pf,ctrVars,d_f);

		for (int j=0; j<(pf+1)*(pf+1); j++)
			ctrVarsFine[i*(pf+1)*(pf+1)+j] = d_f[j];
	}

	// Free variables
	for (int i=0; i<mesh[0].nloc; i++) {
		for (int j=0; j<6; j++)
			delete[] fullBasis[i][j];
		delete[] fullBasis[i];
	}
	delete[] fullBasis;

	for (int i=0; i<(pf+1)*(pf+1); i++) {
		for (int j=0; j<3; j++)
			delete [] fullBasisFine[i][j];
		delete[] fullBasisFine[i];
	}
	delete[] fullBasisFine;

	for (int i=0; i<mesh[0].nloc; i++) {
		for (int j=0; j<6; j++)
			delete [] fullBoundaryBasis[i][j];
		delete[] fullBoundaryBasis[i];
	}
	delete[] fullBoundaryBasis;

	for (int i=0; i<(pf+1)*(pf+1); i++) {
		for (int j=0; j<3; j++)
			delete[] fullBoundaryBasisFine[i][j];
		delete[] fullBoundaryBasisFine[i];
	}
	delete[] fullBoundaryBasisFine;

	// Post-process solution

	// Determine plot points on parent element
	// Assumes same polynomial degree for each element
	// TODO: Add ability to sub-divide parent element further
	int numPlotPts = mesh[0].p[0]+1;
	double interval = 1.0/((double) numPlotPts-1.0);
	double plotPts[numPlotPts];
	for (int i=0; i<numPlotPts; i++)
		plotPts[i] = i*interval;

	// Pre-compute Bernstein basis and derivatives; assumes all elements are of the same type and polynomial degree
	double*** plotBasis;
	plotBasis = new double** [mesh[0].nloc];
	for (int i=0; i<mesh[0].nloc; i++) {
		plotBasis[i] = new double* [3];
		for (int j=0; j<3; j++)
			plotBasis[i][j] = new double [numPlotPts*numPlotPts];
	}
	Full_Bernstein_Basis_and_Derivs(plotBasis,mesh[0].nloc,mesh[0].p[0],numPlotPts,1,plotPts);

	double*** plotBasisFine;
	plotBasisFine = new double** [(pf+1)*(pf+1)];
	for (int i=0; i<(pf+1)*(pf+1); i++) {
		plotBasisFine[i] = new double* [3];
		for (int j=0; j<3; j++)
			plotBasisFine[i][j] = new double [numPlotPts*numPlotPts];
	}
	Full_Bernstein_Basis_and_Derivs(plotBasisFine,(pf+1)*(pf+1),pf,numPlotPts,1,plotPts);

	// Define variables to hold the solution
	double X[numPlotPts*numPlotPts*numElems];
	double Y[numPlotPts*numPlotPts*numElems];
	double U[numPlotPts*numPlotPts*numElems];
	double U_f[numPlotPts*numPlotPts*numElems];

	// Define subelement connectivity array
	int subConnectivity[4*mesh[0].p[0]*mesh[0].p[0]];
	Subelement_Connectivity(subConnectivity,mesh[0].p[0]);

	// Extract the solution
	Extract_Concentration(mesh,plotBasis,numElems,plotPts,numPlotPts,ctrVars,X,Y,U);
	Extract_Fine_Scales(plotBasisFine,numElems,numPlotPts,pf,ctrVarsFine,U_f);

	// Write solution to file
	Write_VTK_File(X,Y,U,U_f,numElems,numPlotPts,subConnectivity,mesh[0].p[0]*mesh[0].p[0]);

	// Free variables
	for (int i=0; i<mesh[0].nloc; i++) {
		for (int j=0; j<3; j++)
			delete[] plotBasis[i][j];
		delete[] plotBasis[i];
	}
	delete[] plotBasis;

	for (int i=0; i<(pf+1)*(pf+1); i++) {
		for (int j=0; j<3; j++)
			delete [] plotBasisFine[i][j];
		delete[] plotBasisFine[i];
	}
	delete[] plotBasisFine;

	PetscFinalize();
	return 0;
}

void Read_IGA_File(ASG_Element*& mesh, int& dim, int& numElems, int& numNodes) {

	string temp, type;
	int elemN, nodeN, nloc, p;
	double pt1, pt2, pt3, wt, extract;
	ifstream inStream;

	// Open IGA file and give error if file cannot be opened
	// (e.g. file does not exist or user does not have write privileges)
	// TODO Edit code so user specifies input file at run time
	inStream.open("test/plate.iga");
	if (inStream.fail())
	{
		cout << "Error: File could not be opened." << endl;
		exit(1);
	}

	// Extract type, number of global nodes, and number of Bezier elements
	inStream >> temp >> type;
	inStream >> temp >> nodeN;
	inStream >> temp >> elemN;

	// Make number of elements and nodes available to main function
	numElems = elemN;
	numNodes = nodeN;

	// Set dimension based on type of geometry
	if (type == "plate") {
		dim = 2;
	} else {
		dim = 3;
	}

	// Initialize arrays to store the global control points and weights
	double points[nodeN*dim];
	double weights[nodeN];

	// Read in and store global control points
	for (int i=0; i<nodeN; i++) {
		inStream >> temp >> pt1 >> pt2 >> pt3 >> wt;

		points[i] = pt1;
		points[nodeN+i] = pt2;
		weights[i] = wt;

		if (dim == 3)
			points[2*nodeN+i] = pt3;
	}

	// Create array of ASG_Elements
	mesh = new ASG_Element[elemN];

	// Loop through elements and extract local control points and weights,
	// connectivity, and extraction operators
	for (int i=0; i<elemN; i++) {

		// Extract number of local basis functions and polynomial order of basis functions
		mesh[i].p = new int[1];
		inStream >> temp >> nloc >> p >> p;
		mesh[i].nloc = nloc;
		mesh[i].p[0] = p;

		// Extract connectivity:
		// Allocate memory for connectivity array
		mesh[i].connectivity = new int[nloc];

		// Read and store connectivity information
		for (int j=0; j<nloc; j++) {
			inStream >> mesh[i].connectivity[j];
		}

		// Extract extraction operator:
		// Initialize array to store extraction operator
		mesh[i].extraction = new double[nloc*nloc];

		// Read and store extraction operator
		for (int j=0; j<nloc; j++) {
			for (int k=0; k<nloc; k++) {
				inStream >> extract;
				mesh[i].extraction[j*nloc+k] = extract;
			}
		}

		// Extract local control points and weights:
		// Initialize arrays to store local control points and weights
		mesh[i].points = new double[nloc*dim];
		mesh[i].weights = new double[nloc];

		// Extract local control points and weights from global control points and weights
		// and store in place holders
		for (int k=0; k<nloc; k++) {
			for (int j=0; j<dim; j++)
				mesh[i].points[k*dim+j] = points[j*nodeN+mesh[i].connectivity[k]];
			mesh[i].weights[k] = weights[mesh[i].connectivity[k]];
		}
	}

	// Close file
	inStream.close();

	return;
}

void Read_BC_File(BC_Struct* boundaryConditions, int* BC, double* g, int numNodes, int numElems) {

	string temp;
	ifstream inStream;

	// Open BC file and give error if file cannot be opened
	// (e.g. file does not exist or user does not have write privileges)
	// TODO Edit code so user specifies input file at run time
	inStream.open("test/plate.bc");
	if (inStream.fail())
	{
		cout << "Error: File could not be opened." << endl;
		exit(1);
	}

	// Extract BC array
	inStream >> temp;
	for (int i=0; i<numNodes; i++) {
		inStream >> BC[i];
	}

	// Extract g array
	inStream >> temp;
	for (int i=0; i<numNodes; i++) {
		inStream >> g[i];
	}

	// Extract Neumann array
	inStream >> temp;
	for (int i=0; i<4; i++) {
		for (int j=0; j<numElems; j++)
			inStream >> boundaryConditions[i].Neumann[j];
	}

	// Extract h array
	inStream >> temp;
	for (int i=0; i<4; i++) {
		for (int j=0; j<numElems; j++)
			inStream >> boundaryConditions[i].h[j];
	}

	// Extract Robin array
	inStream >> temp;
	for (int i=0; i<4; i++) {
		for (int j=0; j<numElems; j++)
			inStream >> boundaryConditions[i].Robin[j];
	}

	// Extract beta array
	inStream >> temp;
	for (int i=0; i<4; i++) {
		for (int j=0; j<numElems; j++)
			inStream >> boundaryConditions[i].beta[j];
	}

	// Close file
	inStream.close();

	return;
}

void Read_Problem_Parameters(Prob_Params* probParams, ASG_Element* mesh, int numNodes, int numElems, double pf) {

	double C[numElems];
	double kappa[numNodes];
	double kappaArt[numNodes];
	double f[numNodes];
	double u[2*numNodes];
	double tempC;
	string temp;
	ifstream inStream;

	// Open VMS file and give error if file cannot be opened
	// (e.g. file does not exist or user does not have write privileges)
	// TODO Edit code so user specifies input file at run time
	inStream.open("test/plate.vms");
	if (inStream.fail())
	{
		cout << "Error: File could not be opened." << endl;
		exit(1);
	}

	// Extract kappa array
	inStream >> temp;
	for (int i=0; i<numNodes; i++) {
		inStream >> kappa[i];
	}

	// Extract kappaArt array
	inStream >> temp;
	for (int i=0; i<numNodes; i++) {
		inStream >> kappaArt[i];
	}

	// Extract f array
	inStream >> temp;
	for (int i=0; i<numNodes; i++) {
		inStream >> f[i];
	}

	// Extract u array
	for (int i=0; i<2; i++) {
		inStream >> temp;
		for (int j=0; j<numNodes; j++) {
			inStream >> u[i*numNodes+j];
		}
	}

	// Extract C array
	inStream >> temp;
	for (int i=0; i<numElems; i++) {
		inStream >> tempC;
		if (tempC == -999)
			C[i] = 8*(pf+1)*(pf+1);
		else
			C[i] = tempC;
	}

	// Loop through elements and store local values for each element
	for (int i=0; i<numElems; i++) {
		for (int j=0; j<mesh[i].nloc; j++) {
			probParams[i].kappa[j] = kappa[mesh[i].connectivity[j]];
			probParams[i].kappaArt[j] = kappaArt[mesh[i].connectivity[j]];
			probParams[i].f[j] = f[mesh[i].connectivity[j]];

			for (int k=0; k<2; k++) {
				probParams[i].u[k*mesh[i].nloc+j] = u[k*numNodes+mesh[i].connectivity[j]];
			}
			probParams[i].C = C[i];
		}
	}

	// Close file
	inStream.close();

	return;
}

void Extract_Geometry(ASG_Element* mesh, int dim, int numElems) {

	int nloc;

	// Create array to hold control points and weights in projected space
	double projPts[nloc*(dim+1)];

	// Create array to hold Bezier control points and weights in projected space
	double bezPts[nloc*(dim+1)];

	// Loop through elements
	for (int i=0; i<numElems; i++) {

		nloc = mesh[i].nloc;

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

		// Initialize arrays to store the Bezier control points and weights
		mesh[i].bezierPoints = new double[nloc*dim];
		mesh[i].bezierWeights = new double[nloc];

		// Assign values to the arrays
		for (int j=0; j<nloc; j++) {
			for (int k=0; k<dim; k++)
				mesh[i].bezierPoints[j*dim+k] = bezPts[j*(dim+1)+k];
			mesh[i].bezierWeights[j] = bezPts[j*(dim+1)+dim];
		}
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
	for (int i=0; i<=n; i++)
		N[i] = new double [p1+1];

	double** M;
	M = new double* [n+1];
	for (int i=0; i<=n; i++)
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
	for (int i=0; i<=n; i++)
		delete[] N[i];
	delete[] N;

	for (int i=0; i<=n; i++)
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
	B = new double* [(p+1)*(p+1)];
	for (int i=0; i<(p+1)*(p+1); i++)
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
	for (int i=0; i<(p+1)*(p+1); i++)
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
	B = new double* [(p+1)*(p+1)];
	for (int i=0; i<(p+1)*(p+1); i++)
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
	for (int i=0; i<(p+1)*(p+1); i++)
		delete[] B[i];
	delete[] B;

	return;
}

void ShapeFunction(ASG_Element element, double** basis, double* R, double** dRdx, double** d2Rdx2, double* x, double* J, int flag) {

	int n_loc = element.nloc;

	for (int i=0; i<n_loc; i++)
		R[i] = 0;
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<2; j++)
			dRdx[i][j] = 0;
	}

	x[0] = 0;
	x[1] = 0;

	double wb;
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
		B[i] = basis[i][0];
	double dBdxi[n_loc][2];
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<2; j++)
			dBdxi[i][j] = basis[i][j+1];
	}

	double d2Bdxi2[n_loc][3];
	if (flag == 1) {
		for (int i=0; i<n_loc; i++) {
			for (int j=0; j<3; j++)
				d2Bdxi2[i][j] = basis[i][j+3];
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

	double d2xdxi2[n_loc][3];
	for (int i=0; i<n_loc; i++) {
		for (int j=0; j<3; j++)
			d2xdxi2[i][j] = 0;
	}

	int indices[n_loc];
	for (int i=0; i<n_loc; i++)
		indices[i] = i;

	double bp[n_loc][2];

	int ptIndices[2];
	ptIndices[0] = 0;
	ptIndices[1] = 1;

	int ind[3];
	for (int j=0; j<3; j++)
		ind[j] = j;

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
				d2Rdxi2[a][1] += element.weights[a]*C_operators[a][b]*((d2Bdxi2[b][1]/wb - dBdxi[b][1]*dwbdxi[1]/(wb*wb)) - (dwbdxi[1]*(dBdxi[b][1]/(wb*wb) - 2*dwbdxi[1]*B[b]/(wb*wb*wb)) + d2wbdxi2[1]*B[b]/(wb*wb)));
				d2Rdxi2[a][2] += element.weights[a]*C_operators[a][b]*((d2Bdxi2[b][2]/wb - dBdxi[b][2]*dwbdxi[2]/(wb*wb)) - (dwbdxi[2]*(dBdxi[b][2]/(wb*wb) - 2*dwbdxi[2]*B[b]/(wb*wb*wb)) + d2wbdxi2[2]*B[b]/(wb*wb)));
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
				d2xdxi2[i][1] += element.bezierWeights[a]*bp[a][i]*((d2Bdxi2[a][1]/wb - dBdxi[a][1]*dwbdxi[1]/(wb*wb)) - (dwbdxi[1]*(dBdxi[a][1]/(wb*wb) - 2*dwbdxi[1]*B[a]/(wb*wb*wb)) + d2wbdxi2[1]*B[a]/(wb*wb)));
				d2xdxi2[i][2] += element.bezierWeights[a]*bp[a][i]*((d2Bdxi2[a][2]/wb - dBdxi[a][2]*dwbdxi[2]/(wb*wb)) - (dwbdxi[2]*(dBdxi[a][2]/(wb*wb) - 2*dwbdxi[2]*B[a]/(wb*wb*wb)) + d2wbdxi2[2]*B[a]/(wb*wb)));
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
				dRdx[a][i] += dRdxi[a][j]*dxidx[j][i];
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

		double b[3];
		double LHS[3];
		int ipiv[3];

		//TODO This could be made more efficient by taking the solve outside the loop?
		for (int i=0; i<n_loc; i++) {
			b[0] = dRdx[i][0]*d2xdxi2[0][0] + dRdx[i][1]*d2xdxi2[1][0];
			b[1] = dRdx[i][0]*d2xdxi2[0][1] + dRdx[i][1]*d2xdxi2[1][1];
			b[2] = dRdx[i][0]*d2xdxi2[0][2] + dRdx[i][1]*d2xdxi2[1][2];

			LHS[0] = d2Rdxi2[i][0] - b[0];
			LHS[1] = d2Rdxi2[i][1] - b[1];
			LHS[2] = d2Rdxi2[i][2] - b[2];

			LAPACKE_dgetrf(LAPACK_ROW_MAJOR,3,3,A,3,ipiv);
			LAPACKE_dgesv(LAPACK_ROW_MAJOR,3,1,A,3,ipiv,LHS,1);

			for (int j=0; j<3; j++)
				d2Rdx2[i][j] = LHS[j];
		}
	}

	return;
}

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

	double** basis;
	basis = new double* [n_loc];
	for (int j=0; j<n_loc; j++)
		basis[j] = new double [6];

	double** basisFine;
	basisFine = new double* [n_locf];
	for (int j=0; j<n_locf; j++)
		basisFine[j] = new double [3];

	double** boundaryBasis;
	boundaryBasis = new double* [n_loc];
	for (int j=0; j<n_loc; j++)
		boundaryBasis[j] = new double [1];

	double** boundaryBasisFine;
	boundaryBasisFine = new double* [n_locf];
	for (int j=0; j<n_locf; j++)
		boundaryBasisFine[j] = new double [3];

	double R[n_loc];

	double** dRdx;
	dRdx = new double* [n_loc];
	for (int j=0; j<n_loc; j++)
		dRdx[j] = new double [2];

	double** d2Rdx2;
	d2Rdx2 = new double* [n_loc];
	for (int j=0; j<n_loc; j++)
		d2Rdx2[j] = new double [3];

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
	h_grid = Grid_Spacing(element.bezierPoints,element.p[0]);

	// Loop through quadrature points

	for (int q1=0; q1<n_q; q1++) {
		for (int q2=0; q2<n_q; q2++) {

			for (int j=0; j<n_loc; j++) {
				for (int k=0; k<6; k++)
					basis[j][k] = fullBasis[j][k][q1*n_q+q2];
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
					diffGalIntegral = k_loc*(dRdx[a][0]*dRdx[b][0] + dRdx[a][1]*dRdx[b][1])*quadWts[q1]*quadWts[q2]*Jdet;
					advGalIntegral = -R[b]*(v_loc[0]*dRdx[a][0] + v_loc[1]*dRdx[a][1])*quadWts[q1]*quadWts[q2]*Jdet;
					k_cc[a*n_loc+b] += diffGalIntegral + advGalIntegral;
				}
				diffTerm = R[a]*f_loc*quadWts[q1]*quadWts[q2]*Jdet;
				f_c[a] += diffTerm;
			}

			for (int af=0; af<n_locf; af++) {
				for (int b=0; b<n_loc; b++){
					advTerm = (v_loc[0]*dRdx[b][0] + v_loc[1]*dRdx[b][1])*Rf[af]*quadWts[q1]*quadWts[q2]*Jdet;
					diffTerm = -k_loc*Rf[af]*(d2Rdx2[b][0] + d2Rdx2[b][2])*quadWts[q1]*quadWts[q2]*Jdet;
					k_fc[af*n_loc+b] += advTerm + diffTerm;
				}
			}

			for (int a=0; a<n_loc; a++) {
				for (int bf=0; bf<n_locf; bf++) {
					k_cf[a*n_locf+bf] += -(v_loc[0]*dRdx[a][0] + v_loc[1]*dRdx[a][1])*Rf[bf]*quadWts[q1]*quadWts[q2]*Jdet;
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
					boundaryBasis[j][k] = fullBoundaryBasis[j][k][side*n_q+q];
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
	for (int j=0; j<n_loc; j++)
		delete [] basis[j];
	delete[] basis;

	for (int j=0; j<n_locf; j++)
		delete[] basisFine[j];
	delete[] basisFine;

	for (int j=0; j<n_loc; j++)
		delete [] boundaryBasis[j];
	delete[] boundaryBasis;

	for (int j=0; j<n_locf; j++)
		delete[] boundaryBasisFine[j];
	delete[] boundaryBasisFine;

	for (int j=0; j<n_loc; j++)
		delete[] dRdx[j];
	delete[] dRdx;

	for (int j=0; j<n_loc; j++)
		delete[] d2Rdx2[j];
	delete[] d2Rdx2;

	for (int j=0; j<n_locf; j++)
		delete[] dRfdx[j];
	delete[] dRfdx;

	return;
}

void Element_Assembly(ASG_Element element, int* BC, double* g, double* k, double* f, Mat K, Vec F) {

	int nloc = element.nloc;
	int i, j;

	for (int a=0; a<nloc; a++){
		i = element.connectivity[a];
		if (BC[i] == 0) {
			for (int b=0; b<nloc; b++) {
				j = element.connectivity[b];
				if (BC[j] == 0) {
					MatSetValue(K,i,j,k[a*nloc+b],ADD_VALUES);
				} else {
					VecSetValue(F,i,-k[a*nloc+b]*g[j],ADD_VALUES);
				}
			}
			VecSetValue(F,i,f[a],ADD_VALUES);
		}
	}

	MatAssemblyBegin(K,MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(K,MAT_FLUSH_ASSEMBLY);

	return;
}

double Get_Local_Value(double vi[], double basis[], int nloc) {

	double locVal = 0;

	for (int i=0; i<nloc; i++)
		locVal += vi[i]*basis[i];

	return locVal;
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

	double** basis;
	basis = new double* [n_loc];
	for (int j=0; j<n_loc; j++)
		basis[j] = new double [6];

	double** basisFine;
	basisFine = new double* [n_locf];
	for (int j=0; j<n_locf; j++)
		basisFine[j] = new double [3];

	double** boundaryBasis;
	boundaryBasis = new double* [n_loc];
	for (int j=0; j<n_loc; j++)
		boundaryBasis[j] = new double [1];

	double** boundaryBasisFine;
	boundaryBasisFine = new double* [n_locf];
	for (int j=0; j<n_locf; j++)
		boundaryBasisFine[j] = new double [3];

	double R[n_loc];

	double** dRdx;
	dRdx = new double* [n_loc];
	for (int j=0; j<n_loc; j++)
		dRdx[j] = new double [2];

	double** d2Rdx2;
	d2Rdx2 = new double* [n_loc];
	for (int j=0; j<n_loc; j++)
		d2Rdx2[j] = new double [3];

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
	h_grid = Grid_Spacing(element.bezierPoints,element.p[0]);

	// Loop through quadrature points

	for (int q1=0; q1<n_q; q1++) {
		for (int q2=0; q2<n_q; q2++) {

			for (int j=0; j<n_loc; j++) {
				for (int k=0; k<6; k++)
					basis[j][k] = fullBasis[j][k][q1*n_q+q2];
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
					advTerm = (v_loc[0]*dRdx[b][0] + v_loc[1]*dRdx[b][1])*Rf[af]*quadWts[q1]*quadWts[q2]*Jdet;
					diffTerm = -k_loc*Rf[af]*(d2Rdx2[b][0] + d2Rdx2[b][2])*quadWts[q1]*quadWts[q2]*Jdet;
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
					boundaryBasis[j][k] = fullBoundaryBasis[j][k][side*n_q+q];
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
	for (int j=0; j<n_loc; j++)
		delete [] basis[j];
	delete[] basis;

	for (int j=0; j<n_locf; j++)
		delete[] basisFine[j];
	delete[] basisFine;

	for (int j=0; j<n_loc; j++)
		delete [] boundaryBasis[j];
	delete[] boundaryBasis;

	for (int j=0; j<n_locf; j++)
		delete[] boundaryBasisFine[j];
	delete[] boundaryBasisFine;

	for (int j=0; j<n_loc; j++)
		delete[] dRdx[j];
	delete[] dRdx;

	for (int j=0; j<n_loc; j++)
		delete[] d2Rdx2[j];
	delete[] d2Rdx2;

	for (int j=0; j<n_locf; j++)
		delete[] dRfdx[j];
	delete[] dRfdx;

	return;
}

void Extract_Concentration(ASG_Element* mesh, double*** fullBasis, int numElems, double* plotPts, int numPlotPts, double* d, double* X, double* Y, double* U) {

	// TODO: Generalize to handle curves in 3D

	// Assume all elements have same nloc
	int nloc = mesh[0].nloc;
	double xi1, xi2, u;

	double** basis;
	basis = new double* [nloc];
	for (int j=0; j<nloc; j++)
		basis[j] = new double [3];

	double R[nloc];

	double** dRdx;
	dRdx = new double* [nloc];
	for (int j=0; j<nloc; j++)
		dRdx[j] = new double [2];

	double** d2Rdx2;
	d2Rdx2 = new double* [nloc];
	for (int j=0; j<nloc; j++)
		d2Rdx2[j] = new double [3];

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
						basis[j][k] = fullBasis[j][k][i1*numPlotPts+i2];
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

	cout << "yo 1889" << endl;
	for (int i=0; i<numPlotPts*numPlotPts*numElems; i++)
		cout << X[i] << " " << Y[i] << " 0" << endl;
	cout << endl;

	// Free variables
	for (int j=0; j<nloc; j++)
		delete [] basis[j];
	delete[] basis;

	for (int j=0; j<nloc; j++)
		delete[] dRdx[j];
	delete[] dRdx;

	for (int j=0; j<nloc; j++)
		delete[] d2Rdx2[j];
	delete[] d2Rdx2;

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

void Subelement_Connectivity(int* subConnectivity, int p) {

	// TODO: Adjust this to incorporate subdivision of subelements
	int numSubElems = p*p;
	int e;

	for (int i=0; i<p; i++) {
		for (int j=0; j<p; j++) {
			e = i*p+j;
			subConnectivity[e] = i*(p+1)+j;
			subConnectivity[e+numSubElems] = i*(p+1)+j+1;
			subConnectivity[e+2*numSubElems] = i*(p+1)+j+(p+1);
			subConnectivity[e+3*numSubElems] = i*(p+1)+j+(p+1)+1;
		}
	}

	return;
}

void Write_VTK_File(double* X, double* Y, double* U, double* U_f, int numElems, int numPlotPts, int* subConnectivity, int numSubElems) {

	ofstream outStream;

	// Open VTK file and give error if file cannot be opened
	// (e.g. directory does not exist or user does not have write privileges)
	// TODO Edit code so user specifies input file at run time
	outStream.open("test/plate.vtk");
	if (outStream.fail())
	{
		cout << "Error: File could not be opened." << endl;
		exit(1);
	}

	// Write file version and identifier
	outStream << "# vtk DataFile Version 3.0\n";

	// TODO Edit code so this is actually useful
	// Write header
	outStream << "Plate example\n";

	// Write file format
	outStream << "ASCII\n\n";

	// Write dataset structure
	outStream << "DATASET UNSTRUCTURED_GRID\n";
	outStream << "POINTS " << numPlotPts*numPlotPts*numElems << " float\n";

	for (int i=0; i<numPlotPts*numPlotPts*numElems; i++)
		outStream << X[i] << " " << Y[i] << " 0\n";
	outStream << "\n";

	outStream << "CELLS " << numElems*numSubElems << " " << (4*numElems*numSubElems) + numElems*numSubElems << "\n";

	for (int i=0; i<numElems; i++) {
		for (int j=0; j<numSubElems; j++) {
			outStream << 4;
			for (int k=0; k<4; k++)
				outStream << " " << i*numPlotPts*numPlotPts+subConnectivity[j*4+k];
			outStream << "\n";
		}
	}
	outStream << "\n";

	outStream << "CELL_TYPES " << numElems*numSubElems << "\n";
	for (int i=0; i<numElems*numSubElems; i++)
		outStream << 8 << "\n";
	outStream << "\n";

	// Write dataset attributes
	outStream << "POINT_DATA " << numPlotPts*numPlotPts*numElems << "\n";
	outStream << "SCALARS concentration float 1\n";
	outStream << "LOOKUP_TABLE default\n";

	for (int i=0; i<numPlotPts*numPlotPts*numElems; i++)
		outStream << U[i] << "\n";

	outStream << "SCALARS fine_scale_sol float 1\n";
	outStream << "LOOKUP_TABLE default\n";

	for (int i=0; i<numPlotPts*numPlotPts*numElems; i++)
		outStream << U_f[i] << "\n";

	// Close file
	outStream.close();

	return;
}
