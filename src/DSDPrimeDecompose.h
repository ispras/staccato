/*******************************************************************************
*   Staccato: Disjoint Support Decompositions from BDDs                        *
*   Copyright (C) 2003-2010  University of Michigan                            *
*   http://www.eecs.umich.edu/staccato                                         *
*                                                                              *
*   Contributors include Stephen Plaza and Valeria Bertacco                    *
*                                                                              *
*   This library is free software; you can redistribute it and/or              *
*   modify it under the terms of the GNU Lesser General Public                 *
*   License as published by the Free Software Foundation; either               *
*   version 2.1 of the License, or (at your option) any later version.         *
*                                                                              *
*   This library is distributed in the hope that it will be useful,            *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of             *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          *
*   Lesser General Public License for more details.                            *
*                                                                              *
*   You should have received a copy of the GNU Lesser General Public           *
*   License along with this library; if not, write to the Free Software        *
*   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA              *
*   02110-1301  USA                                                            *
*                                                                              *
*   Stephen Plaza <splaza@umich.edu>                                           *
*   Valeria Bertacco <valeria@umich.edu>                                       *
*                                                                              *
*   University of Michigan                                                     *
*   Electrical Engineering and Computer Science Dept.                          *
*   2260 Hayward St.                                                           *
*   Ann Arbor, MI 48109-2121                                                   *
*******************************************************************************/


#ifndef DSD_PRIME_DECOMPOSE
#define DSD_PRIME_DECOMPOSE

#include "DSDDecompose.h"

/*!
  Internal function tries different algorithms to see whether the resulting
  DSD will be an inherited prime decomposition
*/
DSDNode* Prime_Decomp(DSDManager* manager, DdNode* f, DdNode *top_func, DSDNode* T, DSDNode* E);

/*!
  Internal function that creates a DSD that is a mux between the variable
  top_func and the functions rooted at the DSDs E and T
*/
DSDNode *BDN_MUX_VAR_DEC_DEC(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *E, DSDNode *T);

/*!
  Internal function that creates a BDD that is a mux between the variable
  top_func and the functions rooted at the DSDs E and T and produces
  the corresponding DSD
*/
DSDNode *BDN_BDD_MUX_VAR_DEC_DEC(DSDManager *manager, DdNode *top_func, DSDNode *E, DSDNode *T);


/*!
  Internal function that creates a BDD whose resulting DSD is
  an inherited prime decomposition
*/
DSDNode *BDN_BDD_PRIME_SUB(DSDManager *manager, DdNode *f, DdNode *top_func, ActualNode *conflict_residue, DSDNode *node);

/*!
  Internal function that creates a BDD whose resulting DSD is
  an inherited prime decomposition
*/
DSDNode *BDN_BDD_PRIME_XOR_SUB(DSDManager *manager, DdNode *f, DdNode *top_func, ActualNode *conflict_residue, DSDNode *node);

/*!
  Internal function that creates a BDD whose resulting DSD is
  an inherited prime decomposition
*/
DSDNode *BDN_BDD_PRIME_MUX_SUB(DSDManager *manager, DdNode *f, DdNode *top_func, ActualNode *conflict_residue, DSDNode *node);

/*!
  Internal function used by the Prime Decomp to identify a potential
  inherited decomposition
*/
ActualNode *cofactor_container_node_equivalence(DSDManager *manager, DSDNode *container, DSDNode *child);

/*!
  Internal function used by the Prime Decomp to identify a potential
  inherited decomposition
*/
ActualNode *cofactor_cofactor_equivalence(DSDManager *manager, DSDNode *container, DSDNode *child, int *switch_them);

/*!
  Internal function used by the Prime Decomp to identify a potential
  inherited decomposition
*/
ActualNode *cofactor_equivalence(DSDManager *manager, DSDNode *container, DSDNode *child);

/*!
  Internal function that creates a prime DSD for a complicated
  special case 
*/
DSDNode *cofactor_elem_BDN_BDD_PRIME_SUB_ELEM(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *node1, DSDNode *node2);

/*!
  Internal function that supports the function cofactor_elem_BDN_BDD_PRIME_SUB_ELEM
*/
int symbolic_finder_builder(DSDManager *manager, DSDNode *node1, DSDNode *node2);

/*!
  Internal function that supports the function cofactor_elem_BDN_BDD_PRIME_SUB_ELEM
*/
int symbolic_finder_builder_incomplete(DSDManager *manager, DSDNode *node1, DSDNode *node2);




#endif
