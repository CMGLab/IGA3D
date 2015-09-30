/*
 * advdiff.h
 *
 *  Created on: Aug 19, 2015
 *      Author: Christopher
 */

#ifndef ADVDIFF_H_
#define ADVDIFF_H_

struct ASG_Element;
struct BC_Struct;
struct Prob_Params;

void Element_Formation(ASG_Element element, BC_Struct elementBC, Prob_Params probParams, double*** fullBasis, double*** fullBasisFine, double*** fullBoundaryBasis, double*** fullBoundaryBasisFine, int n_q, double* quadWts, int pf, double* k_e, double* f_e);

#endif /* ADVDIFF_H_ */
