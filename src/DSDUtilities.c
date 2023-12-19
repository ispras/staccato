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


#include "DSDUtilities.h"

DdNode *symbolic_merger(DdManager *manager, DdNode *base, DdNode *branch, DdNode* top_func)
{
    DdNode *pos, *neg, *result;
   
    assert(base != branch);

    pos = Cudd_Cofactor(manager, base, top_func);
    Cudd_Ref(pos);
    neg = Cudd_Cofactor(manager, base, Cudd_Not(top_func));
    Cudd_Ref(neg);

    result = Cudd_bddIte(manager, branch, pos, neg);	
    Cudd_Ref(result);
    Cudd_RecursiveDeref(manager, pos);
    Cudd_RecursiveDeref(manager, neg);

    Cudd_Deref(result);
    return result;	
}



int_32 __Get_Input_Count(DSDNode *dsd_node)
{
    return INPUT_SIZE(DSD_Regular(dsd_node));
}

DSDNode* __Get_First_Input(DSDNode* dsd_node)
{
    if(dsd_node == NULL)
    {
        return NULL;
    }


    if((DSD_Regular(dsd_node)->actual_list) != NULL)
    {
        return (DSD_Regular(dsd_node)->actual_list)->decomposition;
    }
    else
    {
        return NULL;
    }
}


ActualNode *list_intersection(DdManager *manager, ActualNode *list1, ActualNode *list2, int *size)
{	
    ActualNode *list_iter1;
    ActualNode *list_iter2;
    ActualNode *start;

    ActualNode *final_start;
    ActualNode *temp;

    final_start = NULL;

    start = list2;
    *size = 0;


    for(list_iter1 = list1; list_iter1 != NULL; list_iter1 = list_iter1->next)
    {
        for(list_iter2 = start; list_iter2 != NULL && Cudd_ReadPerm(manager, canonical_var(list_iter2->decomposition)) <= Cudd_ReadPerm(manager, canonical_var(list_iter1->decomposition)); list_iter2 = list_iter2->next)
        {

            start = list_iter2;

            assert(list_iter1->decomposition != NULL);
            assert(list_iter2->decomposition != NULL);

            /*both have to be exactly equal, including sign?*/
            if(list_iter1->decomposition == list_iter2->decomposition)
            {

                start = list_iter2->next;
                (*size)++;

                if(!final_start)
                {
                    final_start = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                    final_start->decomposition = list_iter1->decomposition;
                    final_start->next = NULL;
                    temp = final_start;
                }
                else
                {
                    temp->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                    temp = temp->next;
                    temp->decomposition = list_iter1->decomposition;
                    temp->next = NULL;
                }


                break;	

            }
        }


        if(list_iter2 == NULL)
        {
            break;
        }


    }


    return final_start;


}


void list_intersection_special(DdManager *manager, ActualNode *list1, ActualNode *list2, int *size)
{
    ActualNode *list_iter1;
    ActualNode *list_iter2;
    ActualNode *start;

    ActualNode *final_start;
    ActualNode *temp;

    int opposite_inputs;

    opposite_inputs = 0;


    final_start = NULL;

    start = list2;
    *size = 0;


    for(list_iter1 = list1; list_iter1 != NULL; list_iter1 = list_iter1->next)
    {
        for(list_iter2 = start; list_iter2 != NULL && Cudd_ReadPerm(manager, canonical_var(list_iter2->decomposition)) <= Cudd_ReadPerm(manager, canonical_var(list_iter1->decomposition)); list_iter2 = list_iter2->next)
        {
            start = list_iter2;

            assert(list_iter1->decomposition != NULL);
            assert(list_iter2->decomposition != NULL);

            /*both have to be exactly equal, including sign?*/
            if(list_iter1->decomposition == list_iter2->decomposition)
            {

                start = list_iter2->next;
                (*size)++;			

                break;				
            }
        }

        if(list_iter2 == NULL)
        {
            break;
        }


    }


}



ActualNode *list_residue(DdManager *manager, ActualNode *list1, ActualNode *list2, int *size)
{
    ActualNode *list_iter1;
    ActualNode *list_iter2;

    ActualNode *final_start;
    ActualNode *temp;

    ActualNode *start;

    int mark;


    final_start = NULL;
    mark = 0;


    start = list2;
    (*size) = 0;

    for(list_iter1 = list1; list_iter1 != NULL; list_iter1 = list_iter1->next)
    {
        for(list_iter2 = start; list_iter2 != NULL && Cudd_ReadPerm(manager, canonical_var(list_iter2->decomposition)) <= Cudd_ReadPerm(manager, canonical_var(list_iter1->decomposition)); list_iter2 = list_iter2->next)
        {
            start = list_iter2;

            assert(list_iter1->decomposition != NULL);
            assert(list_iter2->decomposition != NULL);

            /*both have to be exactly equal, including sign?*/
            if(list_iter1->decomposition == list_iter2->decomposition)
            {
                mark = 1;
                start = list_iter2->next;	  
                break;
            }
        }

        if(list_iter2 == NULL)
        {
            start = NULL;
        }

        if(mark == 0)
        {
            (*size)++;

            if(final_start == NULL)
            {
                final_start = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                final_start->decomposition = list_iter1->decomposition;
                final_start->next = NULL;
                temp = final_start;
            }
            else
            {
                temp->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                temp = temp->next;
                temp->decomposition = list_iter1->decomposition;
                temp->next = NULL;
            }
        }

        mark = 0;

    }


    return final_start;

}



int node_exists(ActualNode *list1, DSDNode *node)
{	
    ActualNode *list_iter1;

    assert(node != NULL);

    for(list_iter1 = list1; list_iter1 != NULL; list_iter1 = list_iter1->next)
    {
        if(list_iter1->decomposition == node)
        {
            return 1;				
        }

    }

    return 0;

}


int node_exists_special(ActualNode *list1, DSDNode *node)
{	
    ActualNode *list_iter1;

    assert(node != NULL);

    for(list_iter1 = list1; list_iter1 != NULL; list_iter1 = list_iter1->next)
    {
        if(DSD_Regular(list_iter1->decomposition) == DSD_Regular(node))
        {
            return 1;				
        }

    }

    return 0;

}


/*--------------------------------------------------*/


void set_canonical_var(DSDNode *node)
{
    DSDNode *nodereg;

    ActualNode *iter;

    nodereg = DSD_Regular(node);
    assert(GET_TYPE(nodereg) != VAR);

    iter = nodereg->actual_list;

    while(iter->next != NULL) iter = iter->next;

    SET_CAN((DSD_Regular(nodereg)), (GET_CAN((DSD_Regular(iter->decomposition)))));
}



int canonical_var(DSDNode *node)
{
    DSDNode* nodereg;

    nodereg = DSD_Regular(node);

    return GET_CAN((DSD_Regular(nodereg)));
}


/*rerturn 1 if node1 is greater than node2, 0 if equal,
  -1 if node2 is greater than node1, -2 if disjoint, 2 is some overlap*/ 
int support_compare(DSDManager *manager, DSDNode* node1, DSDNode* node2)
{	
    int size1, size2;
    int max, iter;
    int canonical_var;

    DSDNode *node1_reg;
    DSDNode *node2_reg;

    /*initially disjoint*/

    max = 0;

    node1_reg = DSD_Regular(node1);
    node2_reg = DSD_Regular(node2);

    if(node1_reg->support == NULL)
    {
        support_create(manager, node1_reg);
    }

    if(node2_reg->support == NULL)
    {
        support_create(manager, node2_reg);
    }		

    canonical_var = Cudd_ReadPerm(manager->Ddmanager_analogue, GET_CAN(node1_reg));
    size1 = (canonical_var / (sizeof(int) * 8)) + 1;

    canonical_var = Cudd_ReadPerm(manager->Ddmanager_analogue, GET_CAN(node2_reg)); 
    size2 = (canonical_var / (sizeof(int) * 8)) + 1;

    if(size1 > size2)
    {
        max = size2;
    }
    else
    {
        max = size1;
    }	

    for(iter = 0; iter < max; iter++)
    {
        if((*(node1_reg->support + iter)) & (*(node2_reg->support + iter)))
        {
            return 0;
        }
    }


    return -2;

}

void support_create(DSDManager *manager, DSDNode* node)
{
    SupportList *support;
    SupportList *holder;

    ActualNode *iterActual;

    DSDNode *node_reg;

    node_reg = DSD_Regular(node);

    int canonical_var;
    int support_size;
    int s_iter1;
    int s_iter2;
    int i;

    if(node_reg->support)
    {
        return;
    }

    canonical_var = Cudd_ReadPerm(manager->Ddmanager_analogue, GET_CAN(node_reg));

    support_size = (canonical_var / (sizeof(int) * 8)) + 1;
    manager->support_size += (support_size * sizeof(int));

    if(GET_TYPE(node_reg) == VAR)
    {
        node_reg->support = (int*) malloc(sizeof(int) * support_size);
        memset(node_reg->support, 0, sizeof(int) * support_size);
        (*(node_reg->support + (support_size - 1))) = (*(node_reg->support + ((support_size - 1)))) | (1<<(Cudd_ReadPerm(manager->Ddmanager_analogue, (Cudd_Regular(node_reg->bdd_analogue))->index)%(sizeof(int)*8)));
        return;
    }

    s_iter1 = 0;
    s_iter2 = 0;

    node_reg->support = (int*) malloc(sizeof(int) * support_size);
    memset(node_reg->support, 0, sizeof(int) * support_size);

    iterActual = node_reg->actual_list;

    while(iterActual)
    {


        support_create(manager, iterActual->decomposition);

        canonical_var = Cudd_ReadPerm(manager->Ddmanager_analogue,  GET_CAN((DSD_Regular(iterActual->decomposition))));

        s_iter1 = (canonical_var / (sizeof(int) * 8)) + 1;

        for(i = 0; i < s_iter1; i++)
        {
            (*(node_reg->support + i)) = (*(node_reg->support + i)) | (*(DSD_Regular(iterActual->decomposition)->support + i));
        }

        iterActual = iterActual->next;	
    }


}


DdNode *symbolic_or(DdManager *manager, DSDNode *node)
{
    DdNode *f, *var, *tmp;
    ActualNode *iter;

    assert(!DSD_IsComplement(node));

    f = Cudd_bddIthVar(manager, canonical_var(DSD_Regular(node)->actual_list->decomposition));

    if(DSD_IsComplement(DSD_Regular(node)->actual_list->decomposition))
    {
        f = Cudd_Not(f);
    }

    iter = DSD_Regular(node)->actual_list->next;

    Cudd_Ref(f);


    while(iter)
    {
        var = Cudd_bddIthVar(manager,canonical_var(iter->decomposition));
        
        if(DSD_IsComplement(iter->decomposition))
        {
            var = Cudd_Not(var);
        }

        tmp = Cudd_bddOr(manager,var,f);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager,f);
        f = tmp;

        iter = iter->next;
    }


    cuddDeref(f);
    return f;

}

DdNode *symbolic_xor(DdManager *manager, DSDNode *node)
{
    DdNode *f, *var, *tmp;
    ActualNode *iter;

    f = Cudd_bddIthVar(manager, canonical_var(DSD_Regular(node)->actual_list->decomposition));

    assert(!DSD_IsComplement(node));

    iter = DSD_Regular(node)->actual_list->next;

    Cudd_Ref(f);


    while(iter)
    {
        var = Cudd_bddIthVar(manager,canonical_var(iter->decomposition));

        assert(!DSD_IsComplement(iter->decomposition));

        tmp = Cudd_bddXor(manager,var,f);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager,f);
        f = tmp;

        iter = iter->next;
    }


    cuddDeref(f);
    return f;
}


DdNode *symbolic_mux(DdManager *manager, int top, int e, int t, DSDNode *top_node, DSDNode *Enode, DSDNode *Tnode)
{
    DdNode *f, *gT, *hE;

    f = Cudd_bddIthVar(manager, top);
    if(Tnode && (DSD_Regular(Tnode)->bdd_analogue != Cudd_ReadOne(manager))) gT = Cudd_bddIthVar(manager, t);
    else if(Tnode) gT = Cudd_ReadOne(manager);
    else gT = Cudd_Not(Cudd_ReadOne(manager));

    if(Enode && (DSD_Regular(Enode)->bdd_analogue != Cudd_ReadOne(manager))) hE = Cudd_bddIthVar(manager, e);
    else if(Enode) hE = Cudd_ReadOne(manager);
    else hE = Cudd_Not(Cudd_ReadOne(manager));

    if(DSD_IsComplement(top_node))
    {
        f = Cudd_Not(f);
    }

    if(Enode && DSD_IsComplement(Enode))
    {
        hE = Cudd_Not(hE);
    }

    if(Tnode && DSD_IsComplement(Tnode))
    {
        gT = Cudd_Not(gT);
    }

    f = Cudd_bddIte(manager, f, gT, hE);

    /*Cudd_Ref(f);*/

    return f;
}

void protect(DSDManager *manager, ActualNode *list)
{
    ActualNode *iter;

    iter = list;

    while(iter)
    {
        __DSD_Ref(manager, iter->decomposition);
        iter = iter->next;
    }
}

void unprotect(DSDManager *manager, ActualNode *list)
{
    ActualNode *iter;

    iter = list;

    while(iter)
    {
        __DSD_RecursiveDeref(manager, iter->decomposition);
        iter = iter->next;
    }
}


ActualNode* copy_actual_list(ActualNode *container)
{
    ActualNode *iter1, *iter2, *result;


    iter1 = container;

    result = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    result->decomposition = iter1->decomposition;

    iter2 = result;
    iter1 = iter1->next;


    while(iter1)
    {
        iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

        iter2 = iter2->next;

        iter2->decomposition = iter1->decomposition;
        iter2->next = NULL;

        iter1 = iter1->next;
    }	

    return result;	
}




