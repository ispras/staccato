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


#ifndef DSD_OR_DECOMPOSE
#define DSD_OR_DECOMPOSE

#include "DSDDecompose.h"

/*!
  Internal function tries different algorithms to see whether the resulting
  DSD will be an OR decomposition
*/
DSDNode* OR_Decomp(DSDManager* manager, DdNode* f, DdNode *top_func, DSDNode* T,DSDNode* E);

/*!
  Internal function creates a DSD that is an OR of the variable f and the actuals
  list of base
*/
DSDNode *BDN_OR_VAR_EXP(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base);

/*!
  Internal function creates a DSD that is an NOR of the variable f and the actuals
  list of base
*/
DSDNode *BDN_NOR_VAR_EXP(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base);

/*!
  Internal function creates a DSD that is an OR of the variable f and the function rooted at
  the DSD base
*/
DSDNode *BDN_OR_VAR_DEC(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base);

/*!
  Internal function creates a DSD that is an NOR of the variable f and the function rooted at
  the DSD base
*/
DSDNode *BDN_NOR_VAR_DEC(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base);

/*!
  Internal function creates a DSD that is an OR of the function rooted at node with the remaining
  actuals list being formed from the list actuals
*/
DSDNode *BDN_OR_DEC_ACTUALS(DSDManager *manager, DdNode *f, DSDNode *node, ActualNode *actuals);

/*!
  Internal function creates a DSD that is an NOR of the function rooted at node with the remaining
  actuals list being formed from the list actuals
*/
DSDNode *BDN_NOR_DEC_ACTUALS(DSDManager *manager, DdNode *f, DSDNode *node, ActualNode *actuals);

/*!
  Internal function creates a DSD that is the OR of the functions rooted at node1 and node2 respectively
*/
DSDNode *BDN_OR_DEC_DEC(DSDManager *manager, DdNode *f, DSDNode *node1, DSDNode *node2);

/*!
  Internal function creates a DSD that is the NOR of the functions rooted at node1 and node2 respectively
*/
DSDNode *BDN_NOR_DEC_DEC(DSDManager *manager, DdNode *f, DSDNode *node1, DSDNode *node2);

/*!
  Internal function creates a DSD that is an OR decomposition with the actuals list comprised of actuals
*/
DSDNode *BDN_OR_ACTUALS(DSDManager *manager, DdNode *f, ActualNode *actuals);

/*!
  Internal function creates a DSD that is an NOR decomposition with the actuals list comprised of actuals
*/
DSDNode *BDN_NOR_ACTUALS(DSDManager *manager, DdNode *f, ActualNode *actuals);

/*!
  Internal function creates a BDD of a function that is the OR of the actuals list elements in residue and creates
  a corresponding OR decomposition
*/
DSDNode *BDN_BDD_OR_RESIDUE(DSDManager *manager, ActualNode *residue);

/*!
  Internal function creates a BDD of a function that is the NOR of the actuals list elements in residue and creates
  a corresponding OR decomposition
*/
DSDNode *BDN_BDD_NOR_RESIDUE(DSDManager *manager, ActualNode *residue);

/*!
  Internal function creates a BDD of a function that is the NOR of top_func and the actual list
  elements in residue and creates a corresponding NOR decomposition
*/
DSDNode *BDN_BDD_NOR_VAR_ACTUALS(DSDManager *manager, DdNode *top_func, ActualNode *residue);

/*!
  Internal function creates a BDD of a function that is the NOR of top_func and 
  the function rooted at node and creates a corresponding NOR decomposition
*/
DSDNode *BDN_BDD_NOR_VAR_DEC(DSDManager *manager, DdNode *top_func, DSDNode *node);

#endif
