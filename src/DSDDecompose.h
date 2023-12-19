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


#ifndef DSD_DECOMPOSE
#define DSD_DECOMPOSE

#include "DSD.h"
#include "DSDManager.h"
#include "DSDInterface.h"
#include "DSDOrDecompose.h"
#include "DSDXorDecompose.h"
#include "DSDPrimeDecompose.h"

#include "DSDUtilities.h"

#include "DSDNewDecompose.h"


/*!
  Internal function called by DSD_Create that computes
  the DSD from a BDD, f.  A pointer to the DSDManager
  must be provided. 
*/
DSDNode* __DSD_Create(DSDManager* manager, DdNode* f);

/*!
  Internal function called by __DSD_Create that tries different
  algorithms for identifying the DSD from the BDD.  The variable
  top_func refers to the top BDD variable in f and T and E refer
  to the DSD of the 'then' BDD of f and the 'else' BDD of f
  respectively.
*/
DSDNode* Decomposition(DSDManager* manager, DdNode* f, DdNode *top_func, DSDNode* T, DSDNode* E);

/*!
  Internal function that produces a DSDNode that represents the variable
  corresponding to the BDD variable f.
*/ 
DSDNode *create_var(DSDManager* manager, DdNode* f);

#endif
