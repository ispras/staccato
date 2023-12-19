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


#include "DSDDecompose.h"



void check_one(DSDNode *node)
{
    ActualNode *iter;
    DSDNode *node_reg;
    
    node_reg = DSD_Regular(node);

    if(GET_REF(node_reg) > 0 && GET_TYPE(node_reg) != VAR)
    {
        iter = node_reg->actual_list;
        while(iter)
        {
            check_one(iter->decomposition);
            iter = iter->next;
        }

    }
    else if(GET_TYPE(node_reg) != VAR)
    {
        assert(0);
    }
    
}


/*should caller of DSD_Create reference the new node--why is this done in general.*/

DSDNode* __DSD_Create(DSDManager* manager, DdNode* f)
{
    DSDNode *result;
    DdNode* top_func;
    DdNode* T;
    DdNode* E;

    DSDNode* DT;
    DSDNode* DE;


    DdNode *temp;
    int id;


    if((result = find_DSD_node(manager, f)) != 0)
    {
        return result;
    }

    temp = DD_ONE(manager->Ddmanager_analogue);

    manager->num_entered++;

    /*special case for variables*/
    if((Cudd_Regular(f)->type.kids.T == temp) && (Cudd_Regular(f)->type.kids.E == Cudd_Not(temp)))
    {
        return create_var(manager, f);
    }

    id = Cudd_Regular(f)->index;
    /*temp = manager->Ddmanager_analogue;
      id = *(temp->perm);
      cuddI(manager->Ddmanager_analogue, id);*/	
    top_func = Cudd_ReadVars(manager->Ddmanager_analogue, id);

    if(!Cudd_IsComplement(f))
    {
        T = Cudd_Regular(f)->type.kids.T;
        E = Cudd_Regular(f)->type.kids.E;
    }
    else
    {
        T = Cudd_Not(Cudd_Regular(f)->type.kids.T);
        E = Cudd_Not(Cudd_Regular(f)->type.kids.E);
    }

    DT = __DSD_Create(manager, T);
    
    __DSD_Ref(manager, DT);
    
    DE = __DSD_Create(manager, E);

    __DSD_Ref(manager, DE);
    
    result = Decomposition(manager, f, top_func, DT, DE);

    __DSD_Ref(manager, result);
    
    //assert(check_one(result));
    
    __DSD_RecursiveDeref(manager, DT);
    __DSD_RecursiveDeref(manager, DE);

    //assert(check_one(result));
    
    DSD_Regular(result)->topvar_refsize = ((DSD_Regular(result)->topvar_refsize) & 0xffff0000) | ((DSD_Regular(result)->topvar_refsize & 0x0000ffff) - 1);

    
    /*Decomposition should set type, analogue, actuals, decomp, and next--maybe ref*/

    /*recursively reference?*/

    //check_symbolic2(manager, result);

    return result;

}



DSDNode *create_var(DSDManager* manager, DdNode* f)
{
    DSDNode* result;

    if((result = find_DSD_node(manager, f)) != 0)
    {
        return result;
    }
    
    result = create_DSD_node(manager, Cudd_Regular(f));
    
    //Cudd_RecursiveDeref(manager->Ddmanager_analogue, result->bdd_analogue);
    //manager->num_DSD_nodes--;
    
    assert(Cudd_Regular(f)->index <= SATURATION);
    SET_CAN((DSD_Regular(result)),(Cudd_Regular(f)->index));

    /*DSD node already referenced and initialized)*/
    SET_TYPE(result, VAR);

#ifndef DISABLE_SBDD
    result->symbolic_kernel = Cudd_Regular(f);
    //Cudd_Ref(result->symbolic_kernel);
#endif

    __DSD_Ref(manager, result);
   
    manager->num_DSD_nodes--; 
    if(DSD_IsComplement(f))
    {
        return DSD_Not(result);
    }		


    return result;
}




DSDNode* Decomposition(DSDManager* manager, DdNode* f, DdNode* top_func, DSDNode* T, DSDNode* E)
{
    DSDNode *result;


    if((result = OR_Decomp(manager, f, top_func, T, E)) != 0)
    {
        return result;
    }


    if((result = XOR_Decomp(manager, f, top_func, T, E)) != 0)
    {
        return result;
    }




    if((result = Prime_Decomp(manager, f, top_func, T, E)) != 0)
    {
        manager->num_primes++; 	
        return result;
    } 


    result = Common_Formals_Decomp(manager, f, top_func, T, E);


    return result;

}








