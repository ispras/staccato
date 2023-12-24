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


#include "DSDNewDecompose.h"


static ActualNode *Enode;
static ActualNode *Tnode;

static ActualNode *E_back;
static ActualNode *T_back;


DSDNode* Common_Formals_Decomp(DSDManager* manager, DdNode* f, DdNode *top_func, DSDNode* T, DSDNode* E)
{
    ActualNode *actuals_iter,*actual_holder, *sorted_list, *holder;

    DSDNode *TReg, *EReg;
    DdNode *nodeE_expansion, *nodeT_expansion;
    DSDNode *result;
    int *size;
    int nothing_found;

#ifdef DISABLE_SM
    DdNode *final_expansion, *pos, *neg, *temp, *cube;
    int *length;
#endif

    nothing_found = 0;


    TReg = DSD_Regular(T);
    EReg = DSD_Regular(E);



    /*implement resursive intersection and build*/
    mark_recursive(E, NULL);

    Tnode = NULL;
    Enode = NULL;
    E_back = NULL;
    T_back = NULL;

    /*probably should check to see if T is marked first*/
//    nodeT_expansion = check_marks_recursive(manager, T, NULL);

#ifdef DISABLE_SM    
    if(DSD_IsComplement(T))
    {
        nodeT_expansion = Cudd_Not(TReg->bdd_analogue);
    }
    else
    {
        nodeT_expansion = TReg->bdd_analogue;
    }
#endif

    if(Tnode == NULL)
    {
#ifdef DEBUG_PRINT
        printf("-------Section 4.1----------\n");
#endif

        manager->num_disjoint++;
        unmark_recursive(E);
        return BDN_MUX_VAR_DEC_DEC(manager, f, top_func, E, T);
    }
    else
    {

#ifdef DEBUG_PRINT
        printf("-------Section 4.2----------\n");
#endif

        manager->num_newdecomp += 1;
    }


    //Cudd_Ref(nodeT_expansion);

    nodeE_expansion = check_remaining_marks_recursive(manager, E, NULL);
#ifdef DISABLE_SM
    if(DSD_IsComplement(E))
    {
        nodeE_expansion = Cudd_Not(EReg->bdd_analogue);
    }
    else
    {
        nodeE_expansion = EReg->bdd_analogue;
    }
#endif

    //Cudd_Ref(nodeE_expansion);

    size = (int*) malloc(sizeof(int));
    sorted_list = sort_list(manager->Ddmanager_analogue, size);

#ifdef DISABLE_SM
    #ifndef DISABLE_SBDD
    actuals_iter = sorted_list;
    length = (int*) malloc(sizeof(int));

    while(actuals_iter)
    {
        cube = Cudd_LargestCube(manager->Ddmanager_analogue, actuals_iter->decomposition->bdd_analogue, length);
        if(!cube) continue;
        Cudd_Ref(cube);

        pos = Cudd_Cofactor(manager->Ddmanager_analogue, nodeT_expansion, cube);
        if(!pos) continue;
        Cudd_Ref(pos);

        Cudd_RecursiveDeref(manager->Ddmanager_analogue, cube);
        cube = Cudd_LargestCube(manager->Ddmanager_analogue, Cudd_Not(actuals_iter->decomposition->bdd_analogue), length);
        if(!cube) continue;

        Cudd_Ref(cube);

        neg = Cudd_Cofactor(manager->Ddmanager_analogue, nodeT_expansion, cube);
        if(!neg) continue;
        Cudd_Ref(neg);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, cube);

        temp = Cudd_bddIte(manager->Ddmanager_analogue, Cudd_bddIthVar(manager->Ddmanager_analogue, GET_CAN(DSD_Regular(actuals_iter->decomposition))), pos, neg);	

        Cudd_Ref(temp);

        Cudd_RecursiveDeref(manager->Ddmanager_analogue, nodeT_expansion);	
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, pos);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, neg);

        nodeT_expansion = temp;

        actuals_iter = actuals_iter->next;
    }

    actuals_iter = sorted_list;

    while(actuals_iter)
    {
        cube = Cudd_LargestCube(manager->Ddmanager_analogue, actuals_iter->decomposition->bdd_analogue, length);
        if(!cube) continue;
        Cudd_Ref(cube);


        pos = Cudd_Cofactor(manager->Ddmanager_analogue, nodeE_expansion, cube);
        Cudd_Ref(pos);
        if(!pos) continue;

        Cudd_RecursiveDeref(manager->Ddmanager_analogue, cube);
        cube = Cudd_LargestCube(manager->Ddmanager_analogue, Cudd_Not(actuals_iter->decomposition->bdd_analogue), length);
        if(!cube) continue;
        Cudd_Ref(cube);


        neg = Cudd_Cofactor(manager->Ddmanager_analogue, nodeE_expansion, cube);
        if(!neg) continue;
        Cudd_Ref(neg);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, cube);

        temp = Cudd_bddIte(manager->Ddmanager_analogue, Cudd_bddIthVar(manager->Ddmanager_analogue, GET_CAN(DSD_Regular(actuals_iter->decomposition))), pos, neg);	

        Cudd_Ref(temp);

        Cudd_RecursiveDeref(manager->Ddmanager_analogue, nodeE_expansion);	
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, pos);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, neg);

        nodeE_expansion = temp;

        actuals_iter = actuals_iter->next;
    }


    free(length);

    #endif
#endif

    result = create_DSD_node(manager, f);

    SET_CAN((DSD_Regular(result)), ((Cudd_ReadPerm(manager->Ddmanager_analogue, canonical_var(T)) > Cudd_ReadPerm(manager->Ddmanager_analogue, canonical_var(E))) ? canonical_var(T) : canonical_var(E)));

    assert((*size) > 1);

    (*size)++;	

    SET_TYPE(result, PRIME);
    SET_SIZE(result, *size);
    manager->total_actualsize += *size;

#ifndef DISABLE_SBDD            
    result->symbolic_kernel = Cudd_bddIte(manager->Ddmanager_analogue, top_func, nodeT_expansion, nodeE_expansion);
    Cudd_Ref(result->symbolic_kernel);
    Cudd_RecursiveDeref(manager->Ddmanager_analogue, nodeT_expansion);
    Cudd_RecursiveDeref(manager->Ddmanager_analogue, nodeE_expansion);
#endif
    
    holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);	


    holder->next = sorted_list;
    holder->decomposition = create_var(manager, top_func);

    result->actual_list = holder;


    free(size);

    return result;

}





void insert_listT(DSDManager *manager, ActualNode *actual)
{	
    actual->decomposition = DSD_Regular(actual->decomposition);
    __DSD_Ref(manager, actual->decomposition);

    if(!Tnode)
    {
        actual->next = Tnode;
        Tnode = actual;
        T_back = Tnode;
    }
    else
    {	
        T_back->next = actual;
        actual->next = NULL;
        T_back = T_back->next;
    }
}

void insert_listE(DSDManager *manager, ActualNode *actual)
{	
    actual->decomposition = DSD_Regular(actual->decomposition);
    __DSD_Ref(manager, actual->decomposition);


    if(!Enode)
    {
        actual->next = Enode;
        Enode = actual;
        E_back = Enode;
    }
    else
    {	
        E_back->next = actual;
        actual->next = NULL;
        E_back = E_back->next;
    }
}

void unmark_recursive(DSDNode *node)
{
    ActualNode *actual_iter;

    unmark(DSD_Regular(node));

    DSD_Regular(node)->parent = NULL;

    actual_iter = DSD_Regular(node)->actual_list;

    while(actual_iter)
    {
        unmark_recursive(actual_iter->decomposition);
        actual_iter = actual_iter->next;		
    }

}


ActualNode* sort_list(DdManager *manager, int *size)
{
    ActualNode *result, *iter1, *iter2, *previous, *temp;
    int i, j, marker;

    i = 0;
    j = 0;

    *size = 0;



    result = Tnode;
    iter1 = Tnode;

    if(iter1 == NULL)
    {
        result = Enode;
    }
    else
    {	
        while(iter1->next)
        {
            iter1 = iter1->next;
        }

        iter1->next = Enode;
    }

    if(result == NULL)
    {
        return NULL;
    }

    iter1 = result;

    while(iter1)
    {
        (*size)++;
        iter1 = iter1->next;
    }

    iter1 = result;
    iter2 = result;




    for(i = 0; i < (*size - 1); i++)
    {
        previous = NULL;
        marker = 0;

        for(j = 0, iter2 = result; j < (*size - i - 1); j++)
        {
            if(Cudd_ReadPerm(manager, canonical_var(iter2->decomposition)) > Cudd_ReadPerm(manager, canonical_var(iter2->next->decomposition)))
            {
                marker = 1;
                temp = iter2->next->next;

                if(previous == NULL)
                {
                    result = iter2->next;
                    result->next = iter2;
                    iter2->next = temp;
                    previous = result;
                }
                else
                {
                    previous->next = iter2->next;
                    previous->next->next = iter2;
                    iter2->next = temp;
                    previous = previous->next;
                }
            }
            else
            {
                previous = iter2;
                iter2 = iter2->next;
            }	 

        }		

        if(!marker)
        {
            break;
        }	
    }


    return result;
}


DdNode *check_marks_recursive(DSDManager *manager, DSDNode *node, DSDNode *parent)
{
    ActualNode *actual_iter, *actual_iter2, *or_candidates, *xor_candidates, *temp_candidates, *actual_holder, *simple_candidates, *previous, *actual_insert, *actual_iter3, *actual_iter4, *simple_candidates2;
    DSDNode *parent_node, *dsd_temp;
    DdNode *symbolic_smasher, *node_expansion, *temp;
    int parent_destroyed, elements_remaining, size, parent_type;





    struct cudd_list{
        DdNode *node;
        int top_var;
        struct cudd_list *next;
    };

    struct cudd_list *new_nodes, *temp_node, *node_iter;

    actual_iter = NULL;
    actual_iter2 = NULL;
    actual_iter3 = NULL;
    or_candidates = NULL;
    xor_candidates = NULL;
    actual_holder = NULL;
    simple_candidates = NULL;
    previous = NULL;
    actual_insert = NULL;
    parent_node = NULL;
    symbolic_smasher = NULL;
    dsd_temp = NULL;
    node_expansion = NULL;
    temp = NULL;
    new_nodes = NULL;
    temp_node = NULL;
    node_iter = NULL;
    actual_iter4 = NULL;
    simple_candidates2 = NULL;	

    parent_destroyed = 0;

    actual_iter = DSD_Regular(node)->actual_list;


    if(parent == NULL)
    {
        if(marked(DSD_Regular(node)))
        {
            unmark_recursive(node);

            actual_holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
            actual_holder->decomposition = DSD_Regular(node);
            __DSD_Ref(manager, node);
            actual_holder->next = NULL;
            Tnode = actual_holder;

#ifdef DISABLE_SM
            return NULL;
#else
            node_expansion = Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(node));
            Cudd_Ref(node_expansion);
            if(DSD_IsComplement(node))
            {
                return Cudd_Not(node_expansion);
            }
            else
            {
                return node_expansion;
            }
#endif
        }
    }

    parent_type = GET_TYPE(DSD_Regular(node));


    while(actual_iter)
    {
        if(marked(DSD_Regular(actual_iter->decomposition)))
        {
            parent_destroyed = 1;

            actual_holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
            actual_holder->decomposition = actual_iter->decomposition;



            if((parent_node = DSD_Regular(actual_iter->decomposition)->parent) != NULL)
            {
                if(GET_TYPE(DSD_Regular(parent_node)) == OR && GET_TYPE(DSD_Regular(node)) == OR)
                {
                    if(DSD_IsComplement(parent_node) == DSD_IsComplement(actual_iter->decomposition))
                    {
                        if(or_candidates == NULL)
                        {
                            actual_holder->next = NULL;						 
                            or_candidates = actual_holder;
                            actual_iter2 = or_candidates;
                        }
                        else
                        {
                            actual_iter2->next = actual_holder;
                            actual_iter2 = actual_iter2->next;
                            actual_iter2->next = NULL;
                        }				

                    }
                    else
                    {
                        if(simple_candidates2 == NULL)
                        {
                            actual_holder->next = NULL;						 
                            simple_candidates2 = actual_holder;
                            actual_iter4 = simple_candidates2;
                        }
                        else
                        {
                            actual_iter4->next = actual_holder;
                            actual_iter4 = actual_iter4->next;
                            actual_iter4->next = NULL;
                        }
                    }
                }
                else if(GET_TYPE(DSD_Regular(parent_node)) == XOR && GET_TYPE(DSD_Regular(node)) == XOR)
                {
                    if(xor_candidates == NULL)
                    {
                        actual_holder->next = NULL;						 
                        xor_candidates = actual_holder;
                        actual_iter2 = xor_candidates;
                    }
                    else
                    {
                        actual_iter2->next = actual_holder;
                        actual_iter2 = actual_iter2->next;
                        actual_iter2->next = NULL;
                    }				
                }
                else
                {
                    if(simple_candidates2 == NULL)
                    {
                        actual_holder->next = NULL;						 
                        simple_candidates2 = actual_holder;
                        actual_iter4 = simple_candidates2;
                    }
                    else
                    {
                        actual_iter4->next = actual_holder;
                        actual_iter4 = actual_iter4->next;
                        actual_iter4->next = NULL;
                    }
                }

            }
            else
            {
                if(simple_candidates2 == NULL)
                {
                    actual_holder->next = NULL;						 
                    simple_candidates2 = actual_holder;
                    actual_iter4 = simple_candidates2;
                }
                else
                {
                    actual_iter4->next = actual_holder;
                    actual_iter4 = actual_iter4->next;
                    actual_iter4->next = NULL;
                }
            }		
        }
        else
        {
            if((symbolic_smasher = check_marks_recursive(manager, actual_iter->decomposition, node)) != NULL)
            {
                parent_destroyed = 1;
#ifndef DISABLE_SM
                temp_node = (struct cudd_list *) malloc(sizeof(struct cudd_list));

                temp_node->node = symbolic_smasher;

                SET_CAN((DSD_Regular(temp_node)), (canonical_var(actual_iter->decomposition)));

                temp_node->next = new_nodes;
                new_nodes = temp_node;
#endif
            }
            else
            {

                actual_holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                actual_holder->decomposition = actual_iter->decomposition;

                if(simple_candidates == NULL)
                {
                    actual_holder->next = NULL;						 
                    simple_candidates = actual_holder;
                    actual_iter3 = simple_candidates;
                }
                else
                {
                    actual_iter3->next = actual_holder;
                    actual_iter3 = actual_iter3->next;
                    actual_iter3->next = NULL;
                }
            }
        }
        actual_iter = actual_iter->next;
    }

    if(parent_destroyed)
    {
        actual_iter = simple_candidates;

        elements_remaining = 1;

        if(parent_type == OR)
        {
            node_expansion = Cudd_Not(Cudd_ReadOne(manager->Ddmanager_analogue));
#ifndef DISABLE_SM
            Cudd_Ref(node_expansion);
#endif
            while(elements_remaining)
            {
                size = 0;
                actual_iter = or_candidates;

                if(actual_iter) 
                {
                    parent_node = DSD_Regular(actual_iter->decomposition)->parent;
                }
                else
                {
                    break;
                }

                temp_candidates = actual_iter;

                actual_iter = actual_iter->next;
                or_candidates = actual_iter;

                temp_candidates->next = NULL;
                actual_iter2 = temp_candidates;

                previous = NULL;

                while(actual_iter)
                {
                    if(DSD_Regular(DSD_Regular(actual_iter->decomposition)->parent) == DSD_Regular(parent_node))
                    {
                        actual_iter2->next = actual_iter;
                        actual_iter2 = actual_iter2->next;
                        size++;


                        if(previous == NULL)
                        {
                            actual_iter = actual_iter->next;
                            or_candidates = actual_iter;
                        }
                        else
                        {
                            actual_iter = actual_iter->next;
                            previous->next = actual_iter;
                        }

                        actual_iter2->next = NULL;
                    }
                    else
                    {
                        previous = actual_iter;
                        actual_iter = actual_iter->next;
                    }


                }



                if(size > 0)
                {
                    dsd_temp = BDN_BDD_OR_RESIDUE(manager, temp_candidates);

                    actual_iter = DSD_Regular(dsd_temp)->actual_list;

                    while(actual_iter)
                    {
                        unmark_recursive(actual_iter->decomposition);
                        DSD_Regular(actual_iter->decomposition)->parent = dsd_temp;
                        actual_insert = actual_iter;
                        actual_iter = actual_iter->next;		
                    }

                    temp = node_expansion;

                    actual_holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                    actual_holder->decomposition = dsd_temp;

                    insert_listT(manager, actual_holder);

                    DSD_RecursiveDeref(manager, dsd_temp); 




#ifndef DISABLE_SM
                    node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(dsd_temp)));

                    Cudd_Ref(node_expansion);
                    Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif

                }
                else
                {
                    unmark_recursive(temp_candidates->decomposition);

                    temp = node_expansion;

                    DSD_Regular(temp_candidates->decomposition)->parent = NULL;

#ifndef DISABLE_SM
                    if(DSD_IsComplement(temp_candidates->decomposition))
                    {
                        node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_Not(Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(temp_candidates->decomposition)))); 
                    }
                    else
                    {
                        node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(temp_candidates->decomposition)));
                    }

                    Cudd_Ref(node_expansion);
                    Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif

                    insert_listT(manager, temp_candidates);

                }
            }


            actual_iter = simple_candidates;

            if(actual_iter && actual_iter->next)
            {
                dsd_temp = BDN_BDD_OR_RESIDUE(manager, simple_candidates);

                actual_iter = DSD_Regular(dsd_temp)->actual_list;

                while(actual_iter)
                {
                    actual_insert = actual_iter;
                    actual_iter = actual_iter->next;					
                }

#ifndef DISABLE_SM
                temp = node_expansion;

                node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(dsd_temp)));

                Cudd_Ref(node_expansion);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif

                actual_holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                actual_holder->decomposition = dsd_temp;


                insert_listT(manager, actual_holder);
              
                DSD_RecursiveDeref(manager, dsd_temp); 
 
            
            }
            else
            {
                while(actual_iter)
                {
                    temp = node_expansion;

#ifndef DISABLE_SM
                    if(DSD_IsComplement(actual_iter->decomposition))
                    {
                        node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_Not(Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition)))); 
                    }
                    else
                    {
                        node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition)));
                    }

                    Cudd_Ref(node_expansion);
                    Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif

                    actual_insert = actual_iter;
                    actual_iter = actual_iter->next;


                    insert_listT(manager, actual_insert);
                }
            }		

            actual_iter = simple_candidates2;


            while(actual_iter)
            {
                temp = node_expansion;


                unmark_recursive(actual_iter->decomposition);

                DSD_Regular(actual_iter->decomposition)->parent = NULL;

#ifndef DISABLE_SM
                if(DSD_IsComplement(actual_iter->decomposition))
                {
                    node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_Not(Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition)))); 
                }
                else
                {
                    node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition)));
                }

                Cudd_Ref(node_expansion);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif

                actual_insert = actual_iter;
                actual_iter = actual_iter->next;


                insert_listT(manager, actual_insert);
            }
#ifndef DISABLE_SM

            node_iter = new_nodes;

            while(node_iter)
            {
                temp = node_expansion;

                node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, node_iter->node);

                Cudd_Ref(node_expansion);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, node_iter->node);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);

                temp_node = node_iter;
                node_iter = node_iter->next;
                free(temp_node);
            }

            if(DSD_IsComplement(node))
            {
                node_expansion = Cudd_Not(node_expansion);
            }
#endif
        }
        else if(parent_type == XOR)
        {
            node_expansion = Cudd_Not(Cudd_ReadOne(manager->Ddmanager_analogue));

#ifndef DISABLE_SM
            Cudd_Ref(node_expansion);
#endif
            while(elements_remaining)
            {
                size = 0;
                actual_iter = xor_candidates;

                if(actual_iter) 
                {
                    parent_node = actual_iter->decomposition->parent;
                }
                else
                {
                    break;
                }

                temp_candidates = actual_iter;

                actual_iter = actual_iter->next;
                xor_candidates = actual_iter;

                temp_candidates->next = NULL;
                actual_iter2 = temp_candidates;



                previous = NULL;

                while(actual_iter)
                {
                    assert(!DSD_IsComplement(actual_iter->decomposition));

                    if(actual_iter->decomposition->parent == parent_node)
                    {
                        actual_iter2->next = actual_iter;
                        actual_iter2 = actual_iter2->next;
                        size++;

                        if(previous == NULL)
                        {
                            actual_iter = actual_iter->next;
                            xor_candidates = actual_iter;
                        }
                        else
                        {
                            actual_iter = actual_iter->next;
                            previous->next = actual_iter;
                        }

                        actual_iter2->next = NULL;
                    }
                    else
                    {
                        previous = actual_iter;
                        actual_iter = actual_iter->next;
                    }


                }



                if(size > 0)
                {
                    dsd_temp = BDN_BDD_XOR_RESIDUE(manager, temp_candidates);

                    actual_iter = DSD_Regular(dsd_temp)->actual_list;

                    while(actual_iter)
                    {
                        assert(!DSD_IsComplement(actual_iter->decomposition));
                        unmark_recursive(actual_iter->decomposition);
                        actual_iter->decomposition->parent = dsd_temp;
                        actual_insert = actual_iter;
                        actual_iter = actual_iter->next;										
                    }

#ifndef DISABLE_SM
                    temp = node_expansion;

                    node_expansion = Cudd_bddXor(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(dsd_temp)));

                    Cudd_Ref(node_expansion);
                    Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif
                    actual_holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                    actual_holder->decomposition = dsd_temp;

                    insert_listT(manager, actual_holder);	

                   
                    DSD_RecursiveDeref(manager, dsd_temp); 
 


                }
                else
                {
                    unmark_recursive(temp_candidates->decomposition);

                    insert_listT(manager, temp_candidates);
                    unmark_recursive(temp_candidates->decomposition);

                    DSD_Regular(temp_candidates->decomposition)->parent = NULL;

#ifndef DISABLE_SM                    
                    temp = node_expansion;

                    node_expansion = Cudd_bddXor(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(temp_candidates->decomposition)));

                    Cudd_Ref(node_expansion);
                    Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif

                }

            }


            actual_iter = simple_candidates;
            if(actual_iter && actual_iter->next != NULL)
            {
                dsd_temp = BDN_BDD_XOR_RESIDUE(manager, simple_candidates);

                actual_iter = DSD_Regular(dsd_temp)->actual_list;

                while(actual_iter)
                {
                    actual_insert = actual_iter;
                    actual_iter = actual_iter->next;
                }

#ifndef DISABLE_SM
                temp = node_expansion;

                node_expansion = Cudd_bddXor(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(dsd_temp)));

                Cudd_Ref(node_expansion);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif

                actual_holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                actual_holder->decomposition = dsd_temp;

                insert_listT(manager, actual_holder);
               
                DSD_RecursiveDeref(manager, dsd_temp); 
 
            
            }
            else
            {
                while(actual_iter)
                {

                    temp = node_expansion;

#ifndef DISABLE_SM
                    node_expansion = Cudd_bddXor(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition)));

                    Cudd_Ref(node_expansion);
                    Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif

                    actual_insert = actual_iter;
                    actual_iter = actual_iter->next;
                    insert_listT(manager, actual_insert);

                }
            }

            actual_iter = simple_candidates2;


            while(actual_iter)
            {
                temp = node_expansion;


                unmark_recursive(actual_iter->decomposition);

                DSD_Regular(actual_iter->decomposition)->parent = NULL;

#ifndef DISABLE_SM
                node_expansion = Cudd_bddXor(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition)));

                Cudd_Ref(node_expansion);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif
                actual_insert = actual_iter;
                actual_iter = actual_iter->next;

                insert_listT(manager, actual_insert);


            }
#ifndef DISABLE_SM
            node_iter = new_nodes;

            while(node_iter)
            {
                temp = node_expansion;

                node_expansion = Cudd_bddXor(manager->Ddmanager_analogue, temp, node_iter->node);

                Cudd_Ref(node_expansion);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, node_iter->node);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);

                temp_node = node_iter;
                node_iter = node_iter->next;
                free(temp_node);				
            }


            if(DSD_IsComplement(node))
            {
                node_expansion = Cudd_Not(node_expansion);
            }
#endif
        }
        else
        {
            node_expansion = DSD_Regular(node)->symbolic_kernel;

#ifndef DISABLE_SM
            Cudd_Ref(node_expansion);
#endif
            actual_iter = simple_candidates;

            while(actual_iter)
            {


                actual_insert = actual_iter;
                actual_iter = actual_iter->next;
                insert_listT(manager, actual_insert);

            }

            actual_iter = simple_candidates2;

            while(actual_iter)
            {
                unmark_recursive(actual_iter->decomposition);

                actual_insert = actual_iter;
                actual_iter = actual_iter->next;
                insert_listT(manager, actual_insert);


            }

#ifndef DISABLE_SM
            node_iter = new_nodes;

            while(node_iter)
            {
                temp = node_expansion;

                node_expansion = symbolic_merger(manager->Ddmanager_analogue, temp, node_iter->node, Cudd_bddIthVar(manager->Ddmanager_analogue, GET_CAN((DSD_Regular(node_iter)))));

                Cudd_Ref(node_expansion);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, node_iter->node);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);

                temp_node = node_iter;
                node_iter = node_iter->next;
                free(temp_node);				
            }

            if(DSD_IsComplement(node))
            {
                node_expansion = Cudd_Not(node_expansion);
            }
#endif

        }	
    }
    else
    {
        actual_iter = simple_candidates;		
        while(actual_iter)
        {
            actual_holder = actual_iter;
            actual_iter = actual_iter->next;
            FixHeapFree(actual_malloc_ptr, actual_holder);
        }

        return NULL;
    }	



    return node_expansion;

}



DdNode *check_remaining_marks_recursive(DSDManager *manager, DSDNode *node, DSDNode *parent)
{
    ActualNode *actual_iter, *actual_iter2, *simple_candidates, *actual_holder, *actual_insert, *actual_iter3;
    DSDNode *parent_node, *dsd_temp;
    DdNode *symbolic_smasher, *node_expansion, *temp;
    int parent_destroyed, elements_remaining, size;
    int parent_type, exists;
    int only_one;

    struct cudd_list{
        DdNode *node;
        int top_var;
        struct cudd_list *next;
    };

    struct cudd_list *new_nodes, *temp_node, *node_iter;

    actual_iter = NULL;
    actual_iter2 = NULL;
    actual_iter3 = NULL;
    actual_holder = NULL;
    simple_candidates = NULL;
    actual_insert = NULL;
    parent_node = NULL;
    symbolic_smasher = NULL;
    dsd_temp = NULL;
    node_expansion = NULL;
    temp = NULL;
    new_nodes = NULL;
    temp_node = NULL;
    node_iter = NULL;

    only_one = 1;



    DSD_Regular(node)->parent = NULL;

    parent_destroyed = 0;

    actual_iter = DSD_Regular(node)->actual_list;

    parent_type = GET_TYPE(DSD_Regular(node));


    if(parent == NULL)
    {
        if(!marked(DSD_Regular(node)))
        {
            unmark_recursive(node);		

#ifdef DISABLE_SM
            return 0;
#else
            node_expansion = Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(node));
            Cudd_Ref(node_expansion);
            if(DSD_IsComplement(node))
            {
                return Cudd_Not(node_expansion);
            }
            else
            {
                return node_expansion;
            }
#endif
        }
    }


    unmark(DSD_Regular(node));
    DSD_Regular(node)->parent = NULL;

    while(actual_iter)
    {
        exists = 0;

        if(!marked(DSD_Regular(actual_iter->decomposition)))
        {
            parent_destroyed = 1;
#ifndef DISABLE_SM
            if((parent_node = DSD_Regular(actual_iter->decomposition)->parent) != NULL)
            {
                DSD_Regular(actual_iter->decomposition)->parent = NULL;

                node_iter = new_nodes;

                while(node_iter)
                {
                    if(Cudd_Regular(node_iter->node) == Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(parent_node)))
                    {
                        exists = 1;
                        break;
                    }

                    node_iter = node_iter->next;
                }

                if(!exists)
                {
                    temp_node = (struct cudd_list *) malloc(sizeof(struct cudd_list));
                    temp_node->node = Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(parent_node));
                    Cudd_Ref(temp_node->node);
                    
                    temp_node->next = new_nodes;	
                    new_nodes = temp_node;
                }


            }
            else if(parent_type == OR || parent_type == XOR)
            {
                temp_node = (struct cudd_list *) malloc(sizeof(struct cudd_list));

                if(DSD_IsComplement(actual_iter->decomposition))
                {
                    temp_node->node = Cudd_Not(Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition)));
                    Cudd_Ref(temp_node->node);
                }
                else
                {
                    temp_node->node = Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition));
                    Cudd_Ref(temp_node->node);
                }

                temp_node->next = new_nodes;
                new_nodes = temp_node;
            }		
#endif
        }
        else
        {
            if((symbolic_smasher = check_remaining_marks_recursive(manager, actual_iter->decomposition, node)) != NULL)
            {
                parent_destroyed = 1;
#ifndef DISABLE_SM
                temp_node = (struct cudd_list *) malloc(sizeof(struct cudd_list));
                temp_node->node = symbolic_smasher;


                SET_CAN((DSD_Regular(temp_node)), (canonical_var(actual_iter->decomposition)));

                temp_node->next = new_nodes;
                new_nodes = temp_node;
#endif          
            }
            else
            {
                only_one = 0;

                actual_holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                actual_holder->decomposition = actual_iter->decomposition;

                if(simple_candidates == NULL)
                {
                    actual_holder->next = NULL;						 
                    simple_candidates = actual_holder;
                    actual_iter3 = simple_candidates;
                }
                else
                {
                    actual_iter3->next = actual_holder;
                    actual_iter3 = actual_iter3->next;
                    actual_iter3->next = NULL;
                }

            }
        }


        actual_iter = actual_iter->next;
    }




    if(parent_destroyed)
    {
        actual_iter = simple_candidates;

        elements_remaining = 1;

        if(parent_type == OR)
        {
            node_expansion = Cudd_Not(Cudd_ReadOne(manager->Ddmanager_analogue));
#ifndef DISABLE_SM
            Cudd_Ref(node_expansion);
#endif
            actual_iter = simple_candidates;


            if(actual_iter && actual_iter->next != NULL)
            {
                dsd_temp = BDN_BDD_OR_RESIDUE(manager, simple_candidates);

                actual_iter = DSD_Regular(dsd_temp)->actual_list;

                while(actual_iter)
                {
                    unmark_recursive(actual_iter->decomposition);
                    actual_insert = actual_iter;
                    actual_iter = actual_iter->next;					
                }


#ifndef DISABLE_SM
                temp = node_expansion;


                node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(dsd_temp)));
                Cudd_Ref(node_expansion);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif
                actual_holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                actual_holder->decomposition = dsd_temp;

                insert_listE(manager, actual_holder);					

              
                DSD_RecursiveDeref(manager, dsd_temp); 
 

            }
            else
            {
                while(actual_iter)
                {
#ifndef DISABLE_SM
                    temp = node_expansion;

                    if(DSD_IsComplement(actual_iter->decomposition))
                    {
                        node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_Not(Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition)))); 
                    }
                    else
                    {
                        node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition)));
                    }

                    Cudd_Ref(node_expansion);
                    Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif

                    actual_insert = actual_iter;
                    actual_iter = actual_iter->next;
                    insert_listE(manager, actual_insert);

                }	

            }
#ifndef DISABLE_SM
            node_iter = new_nodes;

            while(node_iter)
            {
                temp = node_expansion;

                node_expansion = Cudd_bddOr(manager->Ddmanager_analogue, temp, node_iter->node);

                Cudd_Ref(node_expansion);

                Cudd_RecursiveDeref(manager->Ddmanager_analogue, node_iter->node);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);

                temp_node = node_iter;
                node_iter = node_iter->next;
                free(temp_node);				
            }

            if(DSD_IsComplement(node))
            {
                node_expansion = Cudd_Not(node_expansion);
            }
#endif
        }
        else if(parent_type == XOR)
        {
            node_expansion = Cudd_Not(Cudd_ReadOne(manager->Ddmanager_analogue));
#ifndef DISABLE_SM
            Cudd_Ref(node_expansion);
#endif
            actual_iter = simple_candidates;

            if(actual_iter && actual_iter->next != NULL)
            {
                dsd_temp = BDN_BDD_XOR_RESIDUE(manager, simple_candidates);

                actual_iter = DSD_Regular(dsd_temp)->actual_list;

                while(actual_iter)
                {
                    actual_insert = actual_iter;
                    actual_iter = actual_iter->next;					
                }
#ifndef DISABLE_SM
                temp = node_expansion;

                node_expansion = Cudd_bddXor(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(dsd_temp)));

                Cudd_Ref(node_expansion);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif
                actual_holder = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
                actual_holder->decomposition = dsd_temp;

                insert_listE(manager, actual_holder);					
              
                DSD_RecursiveDeref(manager, dsd_temp); 
 
                
            }
            else
            {
                while(actual_iter)
                {
                    temp = node_expansion;
#ifndef DISABLE_SM
                    node_expansion = Cudd_bddXor(manager->Ddmanager_analogue, temp, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actual_iter->decomposition)));

                    Cudd_Ref(node_expansion);
                    Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
#endif
                    actual_insert = actual_iter;
                    actual_iter = actual_iter->next;
                    insert_listE(manager, actual_insert);
                }

            }
#ifndef DISABLE_SM
            node_iter = new_nodes;

            while(node_iter)
            {
                temp = node_expansion;

                node_expansion = Cudd_bddXor(manager->Ddmanager_analogue, temp, node_iter->node);

                Cudd_Ref(node_expansion);

                Cudd_RecursiveDeref(manager->Ddmanager_analogue, node_iter->node);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);

                temp_node = node_iter;
                node_iter = node_iter->next;
                free(temp_node);				
            }

            if(DSD_IsComplement(node))
            {
                node_expansion = Cudd_Not(node_expansion);
            }
#endif
        }
        else
        {

            node_expansion = DSD_Regular(node)->symbolic_kernel;

#ifndef DISABLE_SM
            Cudd_Ref(node_expansion);
#endif
            actual_iter = simple_candidates;

            while(actual_iter)
            {
                actual_insert = actual_iter;
                actual_iter = actual_iter->next;
                insert_listE(manager, actual_insert);
            }
#ifndef DISABLE_SM
            node_iter = new_nodes;

            while(node_iter)
            {
                temp = node_expansion;

                node_expansion = symbolic_merger(manager->Ddmanager_analogue, temp, node_iter->node, Cudd_bddIthVar(manager->Ddmanager_analogue, GET_CAN((DSD_Regular(node_iter)))));

                Cudd_Ref(node_expansion);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, node_iter->node);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);

                temp_node = node_iter;
                node_iter = node_iter->next;
                free(temp_node);				
            }

            if(DSD_IsComplement(node))
            {
                node_expansion = Cudd_Not(node_expansion);
            }
#endif
        }	
    }
    else
    {
        actual_iter = simple_candidates;		
        while(actual_iter)
        {
            actual_holder = actual_iter;
            actual_iter = actual_iter->next;
            FixHeapFree(actual_malloc_ptr, actual_holder);
        }


        return NULL;
    }	


    return node_expansion;
}



void mark_recursive(DSDNode *node, DSDNode *parent)
{
    ActualNode *actual_iter;

    DSD_Regular(node)->parent = NULL;

    if(parent != NULL)
    {
        /*should be negative if child is negative?*/
        if(DSD_IsComplement(node) && (GET_TYPE(DSD_Regular(parent)) == OR || GET_TYPE(DSD_Regular(parent)) == XOR))
        {
            assert(GET_TYPE(DSD_Regular(parent)) == OR);

            DSD_Regular(node)->parent = DSD_Complement(parent);
        }
        else if(GET_TYPE(DSD_Regular(parent)) == OR || GET_TYPE(DSD_Regular(parent)) == XOR)
        {
            node->parent = DSD_Regular(parent);
        }
    }

    actual_iter = DSD_Regular(node)->actual_list;

    while(actual_iter)
    {
        mark_recursive(actual_iter->decomposition, node);
        actual_iter = actual_iter->next;		
    }

    mark(DSD_Regular(node));

}

