/*
 * AD2D.cpp
 *
 *  Created on: July 3, 2015
 *  Last modified on: 12 August, 2015
 *      Author: Christopher Coley
 */

#include <petscksp.h>
#include <iostream>
#include "utilities.h"
#include "IGA.h"
#include "VMS.h"
#include "advdiff.h"
#include <time.h>
using namespace std;

int main() {

	clock_t begin, end;
	double time_spent;
	begin = clock();

	//TODO figure out how to do this properly
	PetscInitialize(0,(char***)0,(char*)0,(char*)0);

	int dim, numElems, numNodes, numQ;
	int numDerivs = 2;

	//TODO Include fine scale data stuff somehow
	int pf = 1;

	// Read initialization data from .iga file
	Initialize_Variables(dim,numElems,numNodes);

	// Initialize global variables
	double* globalPoints;
	globalPoints = (double *) malloc(numNodes*dim*sizeof(double));
	double* globalWeights;
	globalWeights = (double *) malloc(numNodes*sizeof(double));
	double** globalElementPoints;
	globalElementPoints = (double **) malloc(numElems*sizeof(double*));
	double** globalElementWeights;
	globalElementWeights = (double **) malloc(numElems*sizeof(double*));
	double** globalBezierPoints;
	globalBezierPoints = (double **) malloc(numElems*sizeof(double*));
	double** globalBezierWeights;
	globalBezierWeights = (double **) malloc(numElems*sizeof(double*));
	double** globalExtraction;
	globalExtraction = (double **) malloc(numElems*sizeof(double*));
	int** globalConnectivity;
	globalConnectivity = (int **) malloc(numElems*sizeof(int*));

	// Create array of ASG_Elements
	ASG_Element* mesh;
	mesh = (ASG_Element *) malloc(numElems*sizeof(ASG_Element));

	// Read in geometry from .iga file
	Read_IGA_File(mesh,globalPoints,globalWeights,globalElementPoints,globalElementWeights,globalBezierPoints,globalBezierWeights,globalExtraction,globalConnectivity);

	end = clock();
	time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
	printf("\nTime for input is %lf\n\n",time_spent);
	begin = clock();

	// Initialize boundary condition variables
	BC_Struct* boundaryConditions;
	boundaryConditions = (BC_Struct *) malloc(numElems*sizeof(BC_Struct));

	int* BC;
	BC = (int *) malloc(numNodes*sizeof(int));
	double* g;
	g = (double *) malloc(numNodes*sizeof(double));
	int** globalElementNeumann;
	globalElementNeumann = (int **) malloc(numElems*sizeof(int*));
	globalElementNeumann[0] = (int *) malloc(numElems*4*sizeof(int));
	double** globalElementh;
	globalElementh = (double **) malloc(numElems*sizeof(double*));
	globalElementh[0] = (double *) malloc(numElems*4*sizeof(double));
	int** globalElementRobin;
	globalElementRobin = (int **) malloc(numElems*sizeof(int*));
	globalElementRobin[0] = (int *) malloc(numElems*4*sizeof(int));
	double** globalElementBeta;
	globalElementBeta = (double **) malloc(numElems*sizeof(double*));
	globalElementBeta[0] = (double *) malloc(numElems*4*sizeof(double));

	for(int i=0; i<numElems; i++) {
		globalElementNeumann[i] = *globalElementNeumann + 4*i;
		globalElementh[i] = *globalElementh + 4*i;
		globalElementRobin[i] = *globalElementRobin + 4*i;
		globalElementBeta[i] = *globalElementBeta + 4*i;
	}

	// Read in boundary condition data from .bc file
	Read_BC_File(boundaryConditions,BC,g,numNodes,numElems,globalElementNeumann,globalElementh,globalElementRobin,globalElementBeta);

	// Assumes all elements have same nloc
	int nloc = mesh[0].nloc;

	// Initialize problem parameter variables
	Prob_Params* probParams;
	probParams = (Prob_Params *) malloc(numElems*sizeof(Prob_Params));

	double** globalElementKappa;
	globalElementKappa = (double **) malloc(numElems*sizeof(double*));
	globalElementKappa[0] = (double *) malloc(numElems*nloc*sizeof(double));
	double** globalElementKappaArt;
	globalElementKappaArt = (double **) malloc(numElems*sizeof(double*));
	globalElementKappaArt[0] = (double *) malloc(numElems*nloc*sizeof(double));
	double** globalElementf;
	globalElementf = (double **) malloc(numElems*sizeof(double*));
	globalElementf[0] = (double *) malloc(numElems*nloc*sizeof(double));
	double** globalElementu;
	globalElementu = (double **) malloc(numElems*sizeof(double*));
	globalElementu[0] = (double *) malloc(2*numElems*nloc*sizeof(double));

	for(int i=0; i<numElems; i++) {
		globalElementKappa[i] = *globalElementKappa + nloc*i;
		globalElementKappaArt[i] = *globalElementKappaArt + nloc*i;
		globalElementf[i] = *globalElementf + nloc*i;
		globalElementu[i] = *globalElementu + 2*nloc*i;
	}

	// Read in problem parameter information
	Read_Problem_Parameters(probParams,mesh,numNodes,numElems,pf,globalElementKappa,globalElementKappaArt,globalElementf,globalElementu);

	end = clock();
	time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
	printf("Time for parameters is %lf\n\n",time_spent);
	begin = clock();

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
	Extract_Geometry(mesh,dim,numElems,globalBezierPoints,globalBezierWeights);

	/*test code to look at geometry assignments
	for (int i=0; i<numElems; i++) {
		cout << "Element " << i+1 << endl;
		cout << "Nloc = " << mesh[i].nloc << endl;
		cout << "p = " << mesh[i].p << endl;
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
	}*/

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
	if (mesh[0].p > pf)
		numQ = mesh[0].p+1;
	else
		numQ = pf+1;

	// Create variables to hold quadrature points and weights
	double quadPts[numQ];
	double quadWts[numQ];

	// Retrieve quadrature points and weights
	Quadrature_Data(numQ,quadPts,quadWts);

	// Pre-compute Bernstein basis and derivatives; assumes all elements are of the same type and polynomial degree
	double*** fullBasis;
	fullBasis = new double** [nloc];
	for (int i=0; i<nloc; i++) {
		fullBasis[i] = new double* [6];
		for (int j=0; j<6; j++)
			fullBasis[i][j] = new double [numQ*numQ];
	}
	Full_Bernstein_Basis_and_Derivs(fullBasis,nloc,mesh[0].p,numQ,numDerivs,quadPts);

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
	fullBoundaryBasis = new double** [nloc];
	for (int i=0; i<nloc; i++) {
		fullBoundaryBasis[i] = new double* [3];
		for (int j=0; j<3; j++)
			fullBoundaryBasis[i][j] = new double [numQ*4];
	}
	Boundary_Bernstein_Basis_and_Derivs(fullBoundaryBasis,nloc,mesh[0].p,numQ,1,quadPts);

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

	MatCreateSeqAIJ(PETSC_COMM_WORLD,numNodes,numNodes,(nloc+mesh[0].p)*(nloc+mesh[0].p),PETSC_NULL,&K);
	MatSetOption(K,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
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

		/* Output for looking at element stiffness matrix and forcing vector
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

	MatAssemblyBegin(K,MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(K,MAT_FLUSH_ASSEMBLY);

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

	end = clock();
	time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
	printf("Time for assembly is %lf\n\n",time_spent);
	begin = clock();

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
	double* ctrVars;
	ctrVars = (double *) malloc(numNodes*sizeof(double));

	int* indices;
	indices = (int *) malloc(numNodes*sizeof(int));
	for (int i=0; i<numNodes; i++)
		indices[i] = i;

	VecGetValues(d,numNodes,indices,ctrVars);

	/* Write out K to a file
	double* Kvals;
	Kvals = new double [numNodes*numNodes];
	MatGetValues(K,numNodes,indices,numNodes,indices,Kvals);
	ofstream outStream;
	outStream.open("test/stiffness.csv");
	for (int i=0; i<numNodes; i++) {
		for (int j=0; j<numNodes; j++)
			outStream << Kvals[i*numNodes+j] << ",";
		outStream << "\n";
	}*/

	/* Write out F to a file
	double* Fvals;
	Fvals = new double [numNodes];
	VecGetValues(F,numNodes,indices,Fvals);
	ofstream outStream2;
	outStream2.open("test/forcing.csv");
	for (int i=0; i<numNodes; i++)
		outStream2 << Fvals[i] << "\n";
	*/

	// Free memory
	MatDestroy(&K);
	VecDestroy(&F);
	VecDestroy(&d);
	KSPDestroy(&ksp);

	// Create variables to hold fine scale solution
	double d_f[(pf+1)*(pf+1)];
	double* ctrVarsFine;
	ctrVarsFine = (double *) malloc(numElems*(pf+1)*(pf+1)*sizeof(double));

	// Solve for the fine scales
	for (int i=0; i<numElems; i++) {
		Solve_Fine_Scales(mesh[i],boundaryConditions[i],probParams[i],fullBasis,fullBasisFine,fullBoundaryBasis,fullBoundaryBasisFine,numQ,quadWts,pf,ctrVars,d_f);

		for (int j=0; j<(pf+1)*(pf+1); j++)
			ctrVarsFine[i*(pf+1)*(pf+1)+j] = d_f[j];
	}

	// Free variables
	for (int i=0; i<nloc; i++) {
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

	for (int i=0; i<nloc; i++) {
		for (int j=0; j<3; j++)
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

	end = clock();
	time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
	printf("Time for solve is %lf\n\n",time_spent);
	begin = clock();

	// Post-process solution

	// Determine plot points on parent element
	// Assumes same polynomial degree for each element
	// TODO: Let user specify refinement
	int numRefine = 2;
	int numIntervals = mesh[0].p*pow(2,numRefine);
	int numPlotPts = numIntervals+1;
	int numSubElems = numIntervals*numIntervals;
	double interval = 1.0/((double) numIntervals);
	double plotPts[numPlotPts];
	for (int i=0; i<numPlotPts; i++)
		plotPts[i] = i*interval;

	// Pre-compute Bernstein basis and derivatives; assumes all elements are of the same type and polynomial degree
	double*** plotBasis;
	plotBasis = new double** [nloc];
	for (int i=0; i<nloc; i++) {
		plotBasis[i] = new double* [3];
		for (int j=0; j<3; j++)
			plotBasis[i][j] = new double [numPlotPts*numPlotPts];
	}
	Full_Bernstein_Basis_and_Derivs(plotBasis,nloc,mesh[0].p,numPlotPts,1,plotPts);

	double*** plotBasisFine;
	plotBasisFine = new double** [(pf+1)*(pf+1)];
	for (int i=0; i<(pf+1)*(pf+1); i++) {
		plotBasisFine[i] = new double* [3];
		for (int j=0; j<3; j++)
			plotBasisFine[i][j] = new double [numPlotPts*numPlotPts];
	}
	Full_Bernstein_Basis_and_Derivs(plotBasisFine,(pf+1)*(pf+1),pf,numPlotPts,1,plotPts);

	// Define variables to hold the solution
	double* X;
	X = (double *) malloc(numPlotPts*numPlotPts*numElems*sizeof(double));
	double* Y;
	Y = (double *) malloc(numPlotPts*numPlotPts*numElems*sizeof(double));
	double* U;
	U = (double *) malloc(numPlotPts*numPlotPts*numElems*sizeof(double));
	double* U_f;
	U_f = (double *) malloc(numPlotPts*numPlotPts*numElems*sizeof(double));

	// Define subelement connectivity array
	int* subConnectivity;
	subConnectivity = (int *) malloc(4*numSubElems*sizeof(int));
	Subelement_Connectivity(subConnectivity,numPlotPts,numSubElems);

	// Extract the solution
	Extract_Concentration(mesh,plotBasis,numElems,plotPts,numPlotPts,ctrVars,X,Y,U);
	Extract_Fine_Scales(plotBasisFine,numElems,numPlotPts,pf,ctrVarsFine,U_f);

	// Write solution to file
	Write_VTK_File(X,Y,U,U_f,numElems,numPlotPts,subConnectivity,numSubElems);

	// Free variables
	for (int i=0; i<nloc; i++) {
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

	end = clock();
	time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
	printf("Time for postprocess is %lf\n\n",time_spent);

	//for (int i=0; i<numElems; i++) {
	//	free(mesh[i].points);
	//	free(mesh[i].weights);
	//	free(mesh[i].bezierPoints);
	//	free(mesh[i].bezierWeights);
	//	free(mesh[i].extraction);
	//	free(mesh[i].connectivity);
	//}
	//free(mesh);

	//for (int i=0; i<numElems; i++) {
	//	delete[] boundaryConditions[i].Neumann;
	//	delete[] boundaryConditions[i].h;
	//	delete[] boundaryConditions[i].Robin;
	//	delete[] boundaryConditions[i].beta;
	//}

	//for (int i=0; i<numElems; i++) {
	//	delete[] probParams[i].kappa;
	//	delete[] probParams[i].kappaArt;
	//	delete[] probParams[i].f;
	//	delete[] probParams[i].u;
	//}

	PetscFinalize();

	cout << "Program completed" << endl;
	return 0;
}
