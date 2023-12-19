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


#ifndef DSD_XOR_DECOMPOSE
#define DSD_XOR_DECOMPOSE

#include "DSDDecompose.h"

/*!
  Internal function tries different algorithms to see whether the resulting
  DSD will be an XOR decomposition
*/
DSDNode* XOR_Decomp(DSDManager* manager, DdNode* f, DdNode *top_func, DSDNode* T, DSDNode* E);

/*!
  Internal function creates a DSD that is an XOR of the variable f and the actuals
  list of base
*/
DSDNode *BDN_XOR_VAR_EXP(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base);

/*!
  Internal function creates a DSD that is a NXOR of the variable f and the actuals
  list of base
*/
DSDNode *BDN_NXOR_VAR_EXP(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base);

/*!
  Internal function creates a DSD that is an XOR of the variable f and the function rooted at
  the DSD base
*/
DSDNode *BDN_XOR_VAR_DEC(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base);

/*!
  Internal function creates a DSD that is a NXOR of the variable f and the function rooted at
  the DSD base
*/
DSDNode *BDN_NXOR_VAR_DEC(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base);

/*!
  Internal function creates a DSD that is an XOR of the function rooted at node with the remaining
  actuals list being formed from the list actuals
*/
DSDNode *BDN_XOR_DEC_ACTUALS(DSDManager *manager, DdNode *f, DSDNode *node, ActualNode *actuals);

/*!
  Internal function creates a DSD that is a NXOR of the function rooted at node with the remaining
  actuals list being formed from the list actuals
*/
DSDNode *BDN_NXOR_DEC_ACTUALS(DSDManager *manager, DdNode *f, DSDNode *node, ActualNode *actuals);


/*!
  Internal function creates a DSD that is the XOR of the functions rooted at node1 and node2 respectively
*/
DSDNode *BDN_XOR_DEC_DEC(DSDManager *manager, DdNode *f, DSDNode *node1, DSDNode *node2);

/*!
  Internal function creates a DSD that is the NXOR of the functions rooted at node1 and node2 respectively
*/
DSDNode *BDN_NXOR_DEC_DEC(DSDManager *manager, DdNode *f, DSDNode *node1, DSDNode *node2);

/*!
  Internal function creates a DSD that is an XOR decomposition with the actuals list comprised of actuals
*/
DSDNode *BDN_XOR_ACTUALS(DSDManager *manager, DdNode *f, ActualNode *actuals);

/*!
  Internal function creates a DSD that is a NXOR decomposition with the actuals list comprised of actuals
*/
DSDNode *BDN_NXOR_ACTUALS(DSDManager *manager, DdNode *f, ActualNode *actuals);

/*!
  Internal function creates a BDD of a function that is the XOR of the actuals list elements in residue and creates
  a corresponding OR decomposition
*/
DSDNode *BDN_BDD_XOR_RESIDUE(DSDManager *manager, ActualNode *residue);

/*!
  Internal function creates a BDD of a function that is the NXOR of the actuals list elements in residue and creates
  a corresponding OR decomposition
*/
DSDNode *BDN_BDD_NXOR_RESIDUE(DSDManager *manager, ActualNode *residue);

#endif
