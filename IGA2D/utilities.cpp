/*
 * utilities.cpp
 *
 *  Created on: Aug 19, 2015
 *      Author: Christopher
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "utilities.h"
#include "IGA.h"
using namespace std;

void Initialize_Variables(int& dim, int& numElems, int& numNodes) {

	string temp, type;
	int elemN, nodeN, nloc, p;
	double pt1, pt2, pt3, wt, connect, extract;
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

	// Read in and throw away global control points
	for (int i=0; i<nodeN; i++)
		inStream >> temp >> pt1 >> pt2 >> pt3 >> wt;

	// Loop through elements and extract number of local basis functions and polynomial order;
	// Throw away local control points and weights, connectivity, and extraction operators
	for (int i=0; i<elemN; i++) {

		// Read and throw away number of local basis functions and polynomial order of basis functions
		inStream >> temp >> nloc >> p >> p;

		// Read and throw away connectivity information
		for (int j=0; j<nloc; j++) {
			inStream >> connect;
		}

		// Read and throw away extraction operator
		for (int j=0; j<nloc; j++) {
			for (int k=0; k<nloc; k++)
				inStream >> extract;
		}
	}

	// Close file
	inStream.close();

	return;
}

void Read_IGA_File(ASG_Element* mesh, double* globalPoints, double* globalWeights, double** globalElementPoints, double** globalElementWeights, double** globalBezierPoints, double** globalBezierWeights, double** globalExtraction, int** globalConnectivity) {

	string temp, type;
	int dim, elemN, nodeN, nloc, p;
	double pt1, pt2, pt3, wt;
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

	// Read in type, number of global nodes, and number of Bezier elements
	inStream >> temp >> type;
	inStream >> temp >> nodeN;
	inStream >> temp >> elemN;

	// Set dimension based on type of geometry
	if (type == "plate") {
		dim = 2;
	} else {
		dim = 3;
	}

	// Read in and store global control points
	for (int i=0; i<nodeN; i++) {
		inStream >> temp >> pt1 >> pt2 >> pt3 >> wt;

		//globalPoints[i] = pt1;
		//globalPoints[nodeN+i] = pt2;
		globalPoints[dim*i] = pt1;
		globalPoints[dim*i+1] = pt2;
		globalWeights[i] = wt;

		if (dim == 3)
			//globalPoints[2*nodeN+i] = pt3;
			globalPoints[dim*i+2] = pt3;
	}

	// Loop through elements and extract local control points and weights,
	// connectivity, and extraction operators
	for (int i=0; i<elemN; i++) {

		// Extract number of local basis functions and polynomial order of basis functions
		inStream >> temp >> nloc >> p >> p;
		mesh[i].nloc = nloc;
		mesh[i].p = p;

		// Initialize variables to store data
		// Assumes nloc is the same for all elements
		if (i == 0) {
			globalElementPoints[0] = (double *) malloc(elemN*nloc*dim*sizeof(double));
			globalElementWeights[0] = (double *) malloc(elemN*nloc*sizeof(double));
			globalBezierPoints[0] = (double *) malloc(elemN*nloc*dim*sizeof(double));
			globalBezierWeights[0] = (double *) malloc(elemN*nloc*sizeof(double));
			globalExtraction[0] = (double *) malloc(elemN*nloc*nloc*sizeof(double));
			globalConnectivity[0] = (int *) malloc(elemN*nloc*sizeof(int));

			for(int j=0; j<elemN; j++) {
				globalElementPoints[j] = *globalElementPoints + nloc*dim*j;
				globalElementWeights[j] = *globalElementWeights + nloc*j;
				globalBezierPoints[j] = *globalBezierPoints + nloc*dim*j;
				globalBezierWeights[j] = *globalBezierWeights + nloc*j;
				globalExtraction[j] = *globalExtraction + nloc*nloc*j;
				globalConnectivity[j] = *globalConnectivity + nloc*j;
			}
		}

		// Extract connectivity:
		// Read and store connectivity information
		for (int j=0; j<nloc; j++) {
			inStream >> globalConnectivity[i][j];
		}

		mesh[i].connectivity = globalConnectivity[i];

		// Extract extraction operator:
		// Read and store extraction operator
		for (int j=0; j<nloc; j++) {
			for (int k=0; k<nloc; k++)
				inStream >> globalExtraction[i][j*nloc+k];
		}
		mesh[i].extraction = globalExtraction[i];

		// Extract local control points and weights from global control points and weights
		// and store in place holders
		for (int k=0; k<nloc; k++) {
			for (int j=0; j<dim; j++)
				globalElementPoints[i][k*dim+j] = globalPoints[mesh[i].connectivity[k]*dim+j];
			globalElementWeights[i][k] = globalWeights[mesh[i].connectivity[k]];
		}
		mesh[i].points = globalElementPoints[i];
		mesh[i].weights = globalElementWeights[i];
	}

	// Close file
	inStream.close();

	return;
}

void Read_BC_File(BC_Struct* boundaryConditions, int* BC, double* g, int numNodes, int numElems, int** globalElementNeumann, double** globalElementh, int** globalElementRobin, double** globalElementBeta) {

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
	for (int i=0; i<numElems; i++) {
		for (int j=0; j<4; j++) {
			inStream >> globalElementNeumann[i][j];
		}
		boundaryConditions[i].Neumann = globalElementNeumann[i];
	}

	// Extract h array
	inStream >> temp;
	for (int i=0; i<numElems; i++) {
		for (int j=0; j<4; j++) {
			inStream >> globalElementh[i][j];
		}
		boundaryConditions[i].h = globalElementh[i];
	}

	// Extract Robin array
	inStream >> temp;
	for (int i=0; i<numElems; i++) {
		for (int j=0; j<4; j++) {
			inStream >> globalElementRobin[i][j];
		}
		boundaryConditions[i].Robin = globalElementRobin[i];
	}

	// Extract beta array
	inStream >> temp;
	for (int i=0; i<numElems; i++) {
		for (int j=0; j<4; j++) {
			inStream >> globalElementBeta[i][j];
		}
		boundaryConditions[i].beta = globalElementBeta[i];
	}

	// Close file
	inStream.close();

	return;
}

void Read_Problem_Parameters(Prob_Params* probParams, ASG_Element* mesh, int numNodes, int numElems, double pf, double** globalElementKappa, double** globalElementKappaArt, double** globalElementf, double** globalElementu) {

	double* C;
	C = (double *) malloc(numElems*sizeof(double));
	double* kappa;
	kappa = (double *) malloc(numNodes*sizeof(double));
	double* kappaArt;
	kappaArt = (double *) malloc(numNodes*sizeof(double));
	double* f;
	f = (double *) malloc(numNodes*sizeof(double));
	double* u;
	u = (double *) malloc(2*numNodes*sizeof(double));
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
			globalElementKappa[i][j] = kappa[mesh[i].connectivity[j]];
			globalElementKappaArt[i][j] = kappa[mesh[i].connectivity[j]];
			globalElementf[i][j] = f[mesh[i].connectivity[j]];

			for (int k=0; k<2; k++)
				globalElementu[i][k*mesh[i].nloc+j] = u[k*numNodes+mesh[i].connectivity[j]];
		}

		probParams[i].kappa = globalElementKappa[i];
		probParams[i].kappaArt = globalElementKappaArt[i];
		probParams[i].f = globalElementf[i];

		for (int k=0; k<2; k++)
			probParams[i].u = globalElementu[i];
		probParams[i].C = C[i];
	}

	// Close file
	inStream.close();

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
		//outStream << X[i] << " " << Y[i] << " " << U[i] << "\n";
	outStream << "\n";

	outStream << "CELLS " << numElems*numSubElems << " " << (4*numElems*numSubElems) + numElems*numSubElems << "\n";

	for (int i=0; i<numElems; i++) {
		for (int j=0; j<numSubElems; j++) {
			outStream << 4;
			for (int k=0; k<4; k++)
				outStream << " " << i*numPlotPts*numPlotPts+subConnectivity[j+k*numSubElems];
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
