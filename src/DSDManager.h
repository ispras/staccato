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


#ifndef DSD_MANAGER
#define DSD_MANAGER

#include "DSD.h"
#include "DSDInterface.h"


/*!
  Internal function used to initialize the various
  caches and structures in STACCATO
*/
DSDManager* __DSD_Init(DdManager* manager, int_32 recommendation_size);

/*!
  Internal function used to reclaim memory and destroy
  outstanding data structures in the given STACCATO
  manager
*/
void __DSD_Quit(DSDManager *manager);

/*!
  Internal function used to reference a DSD node
*/
void __DSD_Ref(DSDManager *manager, DSDNode * dsd);

/*!
  Internal function used to dereference a DSD node
*/
void __DSD_RecursiveDeref(DSDManager *manager, DSDNode * dsd);

/*!
  Internal function to dereference current DSD block
  but does not recursive dereference children blocks
  under any circumstances
*/
void __DSD_Deref(DSDManager *manger, DSDNode *dsd);

/*!
  Internal function called by __DSD_RecursiveDeref
  if the current DSD reference count reaches 0
*/
void recursive_deref(DSDManager *manager, DSDNode *dsd);

/*!
  Internal funnction called by __DSD_Init that creates
  the unique table used to store DSD nodes
*/ 
DSDNode ** Create_DSD_Table(int size);

/*!
  Internal function called by __DSD_Quit to destroy
  unique table used to store DSD nodes
*/
void Destroy_DSD_Table(DSDManager *manager, DSDNode **, int size);

void purge_triggered_stat_update(DSDManager *manager);

/*!
  Internal function called periodically to remove dead
  DSD nodes--garbage cleaning
*/
void Purge_Derefs(DSDManager *manager); 

/*!
  Internal function called by Purge_Derefs to reclaim memory
  used to allocate the actuals list of a deleted DSD node
*/
void delete_actual_list(DSDManager *manager, DSDNode *dsd);
void delete_support(DSDNode *dsd);

/*!
  This function finds a DSD node corresponding to the given
  BDD.  If none exists, 0 is returned.
*/
DSDNode * find_DSD_node(DSDManager *manager, DdNode *bdd); 

/*!
  Internal function to create a new blank DSD node and insert this newly
  created DSD node into the unique table
*/
DSDNode * create_DSD_node(DSDManager *manager, DdNode *bdd); 


/*support for DSD hash table resizing */
/*ASSERT ON MAX VAR SIZE OF 100000--used for a sanity check--should be removed*/
/*not implemented yet*/
/*void adjust_unique_table(DSDManager *manager);*/

#endif
