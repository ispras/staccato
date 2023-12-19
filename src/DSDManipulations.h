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


#ifndef DSD_MANIPULATIONS
#define DSD_MANIPULATIONS

#include "DSD.h"
#include "DSDManager.h"

/*!
  Internal function to print the decomposition type with actuals list
  of a particular DSD node
*/
void  __Decomposition_Print(DSDNode* dsd_node);

/*!
  Internal function to print the full decomposition tree rooted
  at that particular DSD node
*/
void __Recursive_Decomposition_Print(DSDNode* dsd_node);

/*!
  Internal function returns symbolic kernel for a particular
  DSD node
*/
DdNode * __Get_Symbolic_Decomposition(DSDNode* dsd_node);

/*!
  Internal function returns the BDD rooted at a particular
  DSD node
*/
DdNode * __Get_BDD(DSDNode* dsd_node);

#endif
