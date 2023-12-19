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


#include"DSDInterface.h"

int main()
{
  DSDNode *dsd, *dsd_temp;
  DdManager *manager;
  DSDManager *dmanager;
  DdNode *bdd, *temp, *temp2, *symbolic_kernel;
  int num_actuals, type, i, *size;

  /*initialize CUDD manager*/
  manager = Cudd_Init(0,4,CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);

  /*Creating a BDD, a((b XOR c) + d)*/
  temp = Cudd_bddIthVar(manager, 1);
  temp2 = Cudd_bddIthVar(manager, 2);
  temp = Cudd_bddXor(manager, temp, temp2);
  Cudd_Ref(temp);
  temp2 = Cudd_bddIthVar(manager, 3);
  temp2 = Cudd_bddOr(manager, temp, temp2);
  Cudd_Ref(temp2);
  Cudd_RecursiveDeref(manager, temp);
  temp = Cudd_bddIthVar(manager, 0);
  bdd = Cudd_bddAnd(manager, temp, temp2);
  Cudd_Ref(bdd);
  Cudd_RecursiveDeref(manager, temp2);

  /*Initialize DSD manager by choosing a starting cache size*/
  dmanager = DSD_Init(manager, Cudd_ReadMemoryInUse(manager)/2);
  /*Create a DSD from a BDD*/
  dsd = DSD_Create(dmanager, bdd);
  /*Always reference after creation!!!*/
  DSD_Ref(dmanager, dsd);
  
  /*gathering simple decomposition information*/
  symbolic_kernel = get_symbolic_kernel(dsd);
  num_actuals = get_num_actuals(dsd);
  type = Get_Type(dsd);
  Recursive_Decomposition_Print(dsd);

  /*iterating through actuals list--make sure to unmark before calling DSD_Create again*/
  for(i = 0; i < num_actuals; i++)
  {
    dsd_temp = Get_X_Input(dsd, i);
    mark_decomposition(dsd_temp);
  }
  for(i = 0; i < num_actuals; i++)
  {
    dsd_temp = Get_X_Input(dsd, i);
    if(is_marked(dsd_temp))
    {
      unmark_decomposition(dsd_temp);       }
  }

  /*Reclaim DSD and CUDD memory*/
  DSD_Quit(dmanager);
  Cudd_Quit(manager);

  return 0;
}
