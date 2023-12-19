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


#include "DSDManipulations.h"


void  __Decomposition_Print(DSDNode* dsd_node)
{
    ActualNode *iter;

    switch(GET_TYPE(DSD_Regular(dsd_node)))
    {
        case OR:
            if(DSD_IsComplement(dsd_node))
            {
                printf("NOR\n");
            }
            else
            {
                printf("OR\n");
            }		
            break;

        case XOR:
            if(DSD_IsComplement(dsd_node))
            {
                printf("XNOR\n");
            }
            else
            {
                printf("XOR\n");
            }		

            break;

        case PRIME:
            printf("PRIME\n");
            break;


        case VAR:
            printf("VAR\n");	

            break;


    }

    printf("%x\n", dsd_node);

    printf("\tSymbolic Bdd\tCorresponding BDD\n");
    printf("\t%x\t\t%x\n\n", __Get_Symbolic_Decomposition(dsd_node), __Get_BDD(dsd_node));

    iter = DSD_Regular(dsd_node)->actual_list;

    if(iter)
    {
        printf("\t\tDSD Node\tSymbolic Bdd\tCorresponding BDD\n");
    }

    while(iter)
    {
        printf("\t\t%x\t\t%x\t\t%x\n", iter->decomposition, __Get_Symbolic_Decomposition(iter->decomposition), __Get_BDD(iter->decomposition));
        iter = iter->next;
    }


}

void __Recursive_Decomposition_Print(DSDNode* dsd_node)
{
    ActualNode *iter;

    assert(GET_REF(DSD_Regular(dsd_node)) > 0 || (GET_TYPE(DSD_Regular(dsd_node)) == VAR));
    __Decomposition_Print(dsd_node);

    iter = DSD_Regular(dsd_node)->actual_list;

    while(iter)
    {
        __Recursive_Decomposition_Print(iter->decomposition);
        iter = iter->next;
    }

}


DdNode * __Get_Symbolic_Decomposition(DSDNode* dsd_node)
{
    if(DSD_IsComplement(dsd_node))
    {
        return Cudd_Not(DSD_Regular(dsd_node)->symbolic_kernel);
    }
    else
    {
        return dsd_node->symbolic_kernel;
    }	
}

DdNode * __Get_BDD(DSDNode* dsd_node)
{
    if(DSD_IsComplement(dsd_node))
    {
        return Cudd_Not(DSD_Regular(dsd_node)->bdd_analogue);
    }
    else
    {
        return dsd_node->bdd_analogue;
    }	
}
