/*
 * utilities.h
 *
 *  Created on: Aug 19, 2015
 *      Author: Christopher
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

struct ASG_Element;
struct BC_Struct;
struct Prob_Params;

void Initialize_Variables(int& dim, int& numElems, int& numNodes);
void Read_IGA_File(ASG_Element* mesh, double* globalPoints, double* globalWeights, double** globalElementPoints, double** globalElementWeights, double** globalBezierPoints, double** globalBezierWeights, double** globalExtraction, int** globalConnectivity);
void Read_BC_File(BC_Struct* BCs, int* BC, double* g, int numNodes, int numElems, int** globalElementNeumann, double** globalElementh, int** globalElementRobin, double** globalElementBeta);
void Read_Problem_Parameters(Prob_Params* probParams, ASG_Element* mesh, int numNodes, int numElems, double pf,  double** globalElementKappa, double** globalElementKappaArt, double** globalElementf, double** globalElementu);
void Write_VTK_File(double* X, double* Y, double* U, double* U_f, int numElems, int numPlotPts, int* subConnectivity, int numSubElems);

#endif /* UTILITIES_H_ */
