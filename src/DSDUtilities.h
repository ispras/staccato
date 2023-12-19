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


#ifndef DSD_UTILITIES
#define DSD_UTILITIES

#include "DSD.h"
#include "DSDManager.h"
#include "DSDInterface.h"

/*!
  Internal function substitutes the branch BDD for the top_func variable in the base
  function
*/
DdNode *symbolic_merger(DdManager *manager, DdNode *base, DdNode *branch, DdNode *top_func);

/*!
  Internal function that returns the number inputs
  for the DSD dsd_node
*/
int_32 __Get_Input_Count(DSDNode* dsd_node);

/*!
  Internal function that returns the DSD node
  corresponding to the canonically first element
  in the actuals list of DSD dsd_node
*/
DSDNode* __Get_First_Input(DSDNode* dsd_node);

/*!
  Internal function to find common actual lists members between list1 and list2
*/
ActualNode *list_intersection(DdManager *manager, ActualNode *list1, ActualNode *list2, int *size);

/*!
  Internal function returns ActualNodes that are in list1 but
  are not in list2, i.e. list1 - list2
*/
ActualNode *list_residue(DdManager *manager, ActualNode *list1, ActualNode *list2, int *size);

/*!
  Internal function that checks whether the DSDNode node occurs in the actuals list list1
*/
int node_exists(ActualNode *list1, DSDNode *node);

/*!
  Internal function that returns the lowest variable in the decomposition
  tree rooted in the DSD node
*/
int canonical_var(DSDNode *node);

/*!
  Internal function that sets the data structure in DSDNode to explicitly
  state what the lowest variable is in the decomposition tree rooted in this
  DSD
*/
void set_canonical_var(DSDNode *node);

/*!
  Internal function that creates a BDD that is an OR of the representative variables
  of the actuals list of the DSDNode node
*/
DdNode *symbolic_or(DdManager *manager, DSDNode *node);

/*!
  Internal function that creates a BDD that is an XOR of the representative variables
  of the actuals list of the DSDNode node
*/
DdNode *symbolic_xor(DdManager *manager, DSDNode *node);

/*!
  Internal function that creates a BDD that is a mux of the representative variables
  of the DSDNodes top_node, Enode, and Tnode
*/ 
DdNode *symbolic_mux(DdManager *manager, int top, int e, int t, DSDNode *top_node, DSDNode *Enode, DSDNode *Tnode);

/*!
  Internal function that calls __DSD_Ref for each DSDNode in the ActualNode list
*/
void protect(DSDManager *manager, ActualNode *list);

/*!
  Internal function that calls __DSD_RecursiveDeref for each DSDNdoe in the ActualNode list
*/
void unprotect(DSDManager *manager, ActualNode *list);

/*!
  Internal function that produces a copy of the actuals list container
*/
ActualNode* copy_actual_list(ActualNode *container);

/*!
  Internal function to find common actual lists members between list1 and list2 and returns
  the number of actuals list members that are common through the pointer size
*/
void list_intersection_special(DdManager *manager, ActualNode *list1, ActualNode *list2, int *size);


int support_compare(DSDManager *manager, DSDNode* node1, DSDNode* node2);
void support_create(DSDManager *manager, DSDNode* node);

#endif
