/*
 * VMS.h
 *
 *  Created on: Aug 19, 2015
 *      Author: Christopher
 */

#ifndef VMS_H_
#define VMS_H_

struct ASG_Element;
struct BC_Struct;
struct Prob_Params;

void ShapeFine(ASG_Element element, int pf, double** basis, double** dRfdx, double* J);
void Solve_Fine_Scales(ASG_Element element, BC_Struct elementBC, Prob_Params probParams, double*** fullBasis, double*** fullBasisFine, double*** fullBoundaryBasis, double*** fullBoundaryBasisFine, int n_q, double* quadWts, int pf, double* d, double* d_f);
void Extract_Fine_Scales(double*** fullBasisFine, int numElems, int numPlotPts, int pf, double* d, double* U);

#endif /* VMS_H_ */
