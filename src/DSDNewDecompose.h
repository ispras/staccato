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


#ifndef DSD_NEW_DECOMPOSE
#define DSD_NEW_DECOMPOSE

#include "DSDDecompose.h"

/*!
  Internal function that obtains performs the algorithm for computing
  new decompositions of a function from the decompositions of its right and left cofactor
*/
DSDNode* Common_Formals_Decomp(DSDManager* manager, DdNode* f, DdNode *top_func, DSDNode* T, DSDNode* E);

/*!
  Internal function that performs a simple insertion sort to put the actuals list created in the
  new decomposition algorithm in the proper canonical order
*/
ActualNode* sort_list(DdManager *manager, int *size);

/*!
  Internal function called by Common_Formals_Decomp to identify exclusive and common
  blocks by examining the decomposition tree of one of the cofactors
*/
DdNode *check_marks_recursive(DSDManager *manager, DSDNode *node, DSDNode *parent);

/*!
  Internal function called by Common_Formals_Decomp to identify exclusive and common
  blocks by examining the decomposition tree of the other cofactor
*/
DdNode *check_remaining_marks_recursive(DSDManager *manager, DSDNode *node, DSDNode *parent);

/*!
  Internal function to mark all DSD nodes rooted at node
*/
void mark_recursive(DSDNode *node, DSDNode *parent);

/*!
  Internal function to unmark all DSD nodes rooted at node
*/
void unmark_recursive(DSDNode *node);

/*!
  Internal function used to by maintain a list of common and exclusive
  blocks obtained by examining the Then decomposition tree
*/
void insert_listT(DSDManager *manager, ActualNode *actual);

/*!
  Internal function used to by maintain a list of common and exclusive
  blocks obtained by examining the Else decomposition tree
*/
void insert_listE(DSDManager *manager, ActualNode *actual);



#endif
