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


#ifndef DSD_INTERFACE
#define DSD_INTERFACE

#include "DSD.h"
#include "DSDManager.h"
#include "DSDDecompose.h"
#include "DSDManipulations.h"
#include "DSDUtilities.h"

/*!
  A default starting cache value for the DSD unique table.  This value
  should get passed into the DSD_Init function when no prior knowledge is
  known about future size requirements.
*/
#define RECOMMENDATION_DEFAULT 1000000


/*!
 set recommendation_size to some factor of the memory
 used to produce the initial BDDs
*/
DSDManager* DSD_Init(DdManager* manager, int_32 recommendation_size);
void DSD_Quit(DSDManager *manager); /*free resources*/

/*!
  Basic DSD creation function.  It takes a pointer to the 
  DSD Manager along with a BDD
*/
DSDNode *DSD_Create(DSDManager* manager, DdNode* f);

/*!
  This function derefs in the same manner as CUDD does for BDDs
*/
void DSD_RecursiveDeref(DSDManager *manger, DSDNode * dsd_node);

/*!
  This function is the analogue to CUDD version and must be called
  after creating or producing a DSD
*/
void DSD_Ref(DSDManager *manager, DSDNode *dsd_node);

/*!
  This function just derefs the top block in the decomposition.  It behaves
  similarly to the version CUDD uses for BDDs
*/
void DSD_Deref(DSDManager *manager, DSDNode *dsd_node);



/*!
  This function returns the BDD that points to the symbolic
  kernel
*/ 
DdNode * get_symbolic_kernel(DSDNode* dsd_node); 

/*!
  This function returns the BDD that the DSD represents
*/
DdNode * get_bdd(DSDNode* dsd_node);

/*!
  This function returns the type of decomposition, OR, XOR, PRIME, or VAR.
  OR and XOR decompositions are both associative decompositions.  OR, XOR, PRIME,
  and VAR are #defines in STACCATO
*/
int_32 Get_Type(DSDNode* dsd_node); 

/*!
  Get the number of actuals in the actuals list
*/
int_32 get_num_actuals(DSDNode* dsd_node); /*return -1 on error*/

/*!
  This function returns NULL if the variable, index, is less than
  0, is equal to the size of the actuals list, or the dsd_node is
  actually a variable.
*/
DSDNode* Get_X_Input(DSDNode* dsd_node, int index); /*1st actual function*/

/*!
  Used to mark a decomposition.  Make sure to unmark any decompostion
  before computing another decomposition with DSD_Create.
*/
void mark_decomposition(DSDNode *dsd_node);

/*!
  Used to unmark a decomposition
*/
void unmark_decomposition(DSDNode *dsd_node);

/*!
  This function returns 1 when the DSD is marked and 0
  when the DSD is not marked
*/
int is_marked(DSDNode *dsd_node);

/*!
  prints actuals list and decomposition type */
void Decomposition_Print(DSDNode * dsd_node);

/*!
  prints complete decomposition tree
 */
void Recursive_Decomposition_Print(DSDNode * dsd_node);


void load_arrays(DSDManager *manager, DSDNode *result);
void print_stats(DSDManager* manager);
void count_unique(DSDManager* manager, DSDNode *node);
void update_blocks(DSDManager *manager, DSDNode *result);
void validate_deref(DSDManager *manager);
void check_one(DSDNode *node);





#endif
