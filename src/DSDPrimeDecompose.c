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


#include "DSDPrimeDecompose.h"

DSDNode *BDN_MUX_VAR_DEC_DEC(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *E, DSDNode *T)
{
    DSDNode *top_node, *efake, *tfake;
    DSDNode *result;



    ActualNode *iter2;

    DSDNode *temp;
    int var1, var2, temp_var, evar, tvar, efakevar, tfakevar;


    result = create_DSD_node(manager, f);




    SET_TYPE(result, DSD_PRIME);


    if(E && E != manager->one)	evar = canonical_var(E);
    else evar = -1;

    if(T && T != manager->one) tvar = canonical_var(T);
    else tvar = -1;

    efakevar = evar;
    tfakevar = tvar;

    efake = E;
    tfake = T;

    if(Cudd_ReadPerm(manager->Ddmanager_analogue, evar) > Cudd_ReadPerm(manager->Ddmanager_analogue, tvar))
    {
        temp = E;
        E = T;
        T = temp;


        var1 = tvar;
        var2 = evar;		
    }
    else
    {
        var2 = tvar;
        var1 = evar;
    }

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);
    assert(var1 < 10000);
    assert(var2 < 10000);


    top_node = create_var(manager, top_func);

    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    result->actual_list->decomposition = DSD_Regular(top_node);

    result->actual_list->next = NULL;

    iter2 = result->actual_list;


    if(var1 != -1)
    {
        iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
        iter2 = iter2->next;	

        iter2->decomposition = DSD_Regular(E);


        __DSD_Ref(manager, iter2->decomposition);


    }


    if(var2 != -1)
    {
        iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
        iter2 = iter2->next;

        iter2->decomposition = DSD_Regular(T);
        iter2->next = NULL;

        __DSD_Ref(manager, iter2->decomposition);
    }


#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_mux(manager->Ddmanager_analogue, Cudd_Regular(top_func)->index, efakevar, tfakevar, top_node, efake, tfake);
    Cudd_Ref(result->symbolic_kernel);
#endif
    
    set_canonical_var(result);

    SET_SIZE(result, 3);
    manager->total_actualsize += 3;


    return result;
}




DSDNode *BDN_BDD_MUX_VAR_DEC_DEC(DSDManager *manager, DdNode *top_func, DSDNode *E, DSDNode *T)
{
    DSDNode *result;

    DdNode *f;
    DdNode *cuddE;
    DdNode *cuddT;


    if(E && DSD_IsComplement(E))
    {
        cuddE = Cudd_Not(DSD_Regular(E)->bdd_analogue);
    }
    else if(E)
    {
        cuddE = E->bdd_analogue;
    }
    else
    {
        cuddE = Cudd_Not(Cudd_ReadOne(manager->Ddmanager_analogue));
    }

    if(T && DSD_IsComplement(T))
    {
        cuddT = Cudd_Not(DSD_Regular(T)->bdd_analogue);
    }
    else if(T)
    {
        cuddT = T->bdd_analogue;
    }
    else
    {
        cuddT = Cudd_Not(Cudd_ReadOne(manager->Ddmanager_analogue));
    }

    f = Cudd_bddIte(manager->Ddmanager_analogue, top_func, cuddT, cuddE);

    Cudd_Ref(f);

    /*important change*/	

    result = __DSD_Create(manager, f);
    __DSD_Ref(manager, result);
    
    /*	if((result = find_DSD_node(manager, f)) == NULL)
        {
        result = BDN_MUX_VAR_DEC_DEC(manager, f, top_func, E, T);
        }
     */	

    Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);

    return result;	

}


/*bdd creation for sub value*/
DSDNode *BDN_BDD_PRIME_SUB(DSDManager *manager, DdNode *f, DdNode *top_func, ActualNode *conflict_residue, DSDNode *node)
{
    int negative_cofactor, complemented, complemented_var;

    ActualNode *iter;

    DSDNode *sub, *result, *node_reg;
    DdNode *temp, *func;


    negative_cofactor = 0;

    node_reg = DSD_Regular(node);

    complemented = 0;
    complemented_var = 0;


    if(DSD_IsComplement(conflict_residue))
    {
        negative_cofactor = 1;
        conflict_residue = (ActualNode *) DSD_Regular(conflict_residue);
    }

    iter = conflict_residue;

    while(iter != NULL)
    {

        if(marked(DSD_Regular(iter->decomposition)))
        {
            unmark(DSD_Regular(iter->decomposition));

            assert(!DSD_IsComplement(iter->decomposition));

            /*unique part*/

            if(negative_cofactor)
            {
                temp = iter->decomposition->bdd_analogue;

                func = Cudd_bddAnd(manager->Ddmanager_analogue, top_func, temp);

                Cudd_Ref(func);

                if((sub = find_DSD_node(manager, func)) == NULL)
                {
                    /*if(GET_TYPE(iter->decomposition) != DSD_OR)
                      {*/
                    sub = BDN_NOR_VAR_DEC(manager, func, Cudd_Not(top_func), DSD_Not(iter->decomposition));
                    /*}
                      else
                      {
                      sub = BDN_NOR_VAR_EXP(manager, func, Cudd_Not(top_func), DSD_Not(iter->decomposition));
                      }*/

                }


            }
            else
            {
                temp = iter->decomposition->bdd_analogue;

                func = Cudd_bddOr(manager->Ddmanager_analogue, Cudd_Not(top_func), temp);

                Cudd_Ref(func);

                if((sub = find_DSD_node(manager, func)) == NULL)
                {
                    if(GET_TYPE(iter->decomposition) != DSD_OR)
                    {
                        sub = BDN_OR_VAR_DEC(manager, func, Cudd_Not(top_func), iter->decomposition);
                    }
                    else
                    {
                        sub = BDN_OR_VAR_EXP(manager, func, Cudd_Not(top_func), iter->decomposition);
                    }


                }

            }

            Cudd_RecursiveDeref(manager->Ddmanager_analogue, func);



            /*end unique part*/

            iter->decomposition = DSD_Regular(sub);


            if(DSD_IsComplement(sub))
            {
                complemented = 1;
                complemented_var = canonical_var(iter->decomposition);
            }

        }



        DSD_Ref(manager, iter->decomposition);



        iter = iter->next;
    }

    result = create_DSD_node(manager, f);

    __DSD_Ref(manager, result);
    
    result->type_actualsize = node_reg->type_actualsize;
    
    
    manager->total_actualsize += INPUT_SIZE(DSD_Regular(result));
    
    SET_CAN((DSD_Regular(result)), (canonical_var(node_reg)));

    result->actual_list = conflict_residue;


#ifndef DISABLE_SBDD
    result->symbolic_kernel = node_reg->symbolic_kernel;
    Cudd_Ref(result->symbolic_kernel);
    
    if(complemented)
    {
        temp = result->symbolic_kernel;
        result->symbolic_kernel = symbolic_merger(manager->Ddmanager_analogue, temp, Cudd_Not(Cudd_bddIthVar(manager->Ddmanager_analogue, complemented_var)), Cudd_bddIthVar(manager->Ddmanager_analogue, complemented_var));
        Cudd_Ref(result->symbolic_kernel);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
    }
#endif
    

    if(DSD_IsComplement(node))
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        result = DSD_Not(result);
    }

    return result;

}


DSDNode *BDN_BDD_PRIME_XOR_SUB(DSDManager *manager, DdNode *f, DdNode *top_func, ActualNode *conflict_residue, DSDNode *node)
{
    int negative_cofactor;
    int complemented, complemented_var;


    ActualNode *iter;

    DSDNode *sub, *result, *node_reg;

    DdNode *temp, *func;

    negative_cofactor = 0;

    complemented = 0;
    complemented_var = 0;

    node_reg = DSD_Regular(node);



    iter = conflict_residue;

    while(iter != NULL)
    {

        if(marked(DSD_Regular(iter->decomposition)))
        {
            unmark(DSD_Regular(iter->decomposition));

            /*unique part*/


            if(DSD_IsComplement(iter->decomposition))
            {
                temp = Cudd_Not(DSD_Regular(iter->decomposition)->bdd_analogue);
            }
            else
            {
                temp = iter->decomposition->bdd_analogue;
            }

            func = Cudd_bddXor(manager->Ddmanager_analogue, top_func, temp);

            func = Cudd_Not(func);

            Cudd_Ref(func);

            if((sub = find_DSD_node(manager, func)) == NULL)
            {
                if(GET_TYPE(iter->decomposition) != DSD_XOR)
                {
                    sub = BDN_NXOR_VAR_DEC(manager, func, top_func, iter->decomposition);
                }
                else
                {
                    sub = BDN_NXOR_VAR_EXP(manager, func, top_func, iter->decomposition);
                }
            }


            Cudd_RecursiveDeref(manager->Ddmanager_analogue, func);



            /*end unique part*/

            iter->decomposition = DSD_Regular(sub);

            if(DSD_IsComplement(sub))
            {
                complemented = 1;
                complemented_var = canonical_var(iter->decomposition);
            }


        }



        DSD_Ref(manager, iter->decomposition);



        iter = iter->next;
    }

    result = create_DSD_node(manager, f);

    __DSD_Ref(manager, result);
    
    result->type_actualsize = node_reg->type_actualsize;

    manager->total_actualsize += INPUT_SIZE(DSD_Regular(result));
    
    SET_CAN((DSD_Regular(result)), (canonical_var(node_reg)));

    result->actual_list = conflict_residue;

#ifndef DISABLE_SBDD
    result->symbolic_kernel = node_reg->symbolic_kernel;
    Cudd_Ref(result->symbolic_kernel);

    if(complemented)
    {
        temp = result->symbolic_kernel;
        result->symbolic_kernel = symbolic_merger(manager->Ddmanager_analogue, temp, Cudd_Not(Cudd_bddIthVar(manager->Ddmanager_analogue, complemented_var)), Cudd_bddIthVar(manager->Ddmanager_analogue, complemented_var));
        Cudd_Ref(result->symbolic_kernel);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
    }
#endif
    
    if(DSD_IsComplement(node))
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        result = DSD_Not(result);
    }

    return result;
}



/*MODIFY TO ALLOW EXTRA CALCULATION OF FIRST ELEMENT WHICH IS NECESSARILY THE Z*/
DSDNode *BDN_BDD_PRIME_MUX_SUB(DSDManager *manager, DdNode *f, DdNode *top_func, ActualNode *conflict_residue, DSDNode *node)
{
    int negative_cofactor, complemented, complemented_var;

    ActualNode *iter;

    DSDNode *sub, *result, *node_reg;
    DSDNode *Tnode;

    DdNode *temp, *tempT, *func;

    negative_cofactor = 0;

    node_reg = DSD_Regular(node);


    complemented = 0;
    complemented_var = 0;

    


    if(DSD_IsComplement(conflict_residue))
    {
        negative_cofactor = 1;
        conflict_residue = (ActualNode *) DSD_Regular(conflict_residue);

        Tnode = DSD_Not(conflict_residue->decomposition);	
    }
    else
    {
        Tnode = conflict_residue->decomposition;
    }


    protect(manager, conflict_residue->next);




    iter = conflict_residue->next;

    while(iter != NULL)
    {

        if(marked(DSD_Regular(iter->decomposition)))
        {
            unmark(DSD_Regular(iter->decomposition));


            if(negative_cofactor)
            {
                if(DSD_IsComplement(iter->decomposition))
                {
                    temp = Cudd_Not(DSD_Regular(iter->decomposition)->bdd_analogue);
                }
                else
                {
                    temp = iter->decomposition->bdd_analogue;
                }


                if(DSD_IsComplement(Tnode))
                {
                    tempT = Cudd_Not(DSD_Regular(Tnode)->bdd_analogue);
                }
                else
                {
                    tempT = Tnode->bdd_analogue;
                }


                func = Cudd_bddIte(manager->Ddmanager_analogue, top_func, tempT, temp);

                Cudd_Ref(func);

                if((sub = find_DSD_node(manager, func)) == NULL)
                {
                    sub = __DSD_Create(manager, func);
                }


            }
            else
            {

                if(DSD_IsComplement(iter->decomposition))
                {
                    temp = Cudd_Not(DSD_Regular(iter->decomposition)->bdd_analogue);
                }
                else
                {
                    temp = iter->decomposition->bdd_analogue;
                }

                if(DSD_IsComplement(Tnode))
                {
                    tempT = Cudd_Not(DSD_Regular(Tnode)->bdd_analogue);
                }
                else
                {
                    tempT = Tnode->bdd_analogue;
                }

                func = Cudd_bddIte(manager->Ddmanager_analogue, top_func, tempT, temp);

                Cudd_Ref(func);

                if((sub = find_DSD_node(manager, func)) == NULL)
                {
                    sub = __DSD_Create(manager, func);
                    /*sub = BDN_MUX_VAR_DEC_DEC(manager, func, top_func, iter->decomposition, Tnode);*/
                }

            }

            Cudd_RecursiveDeref(manager->Ddmanager_analogue, func);

            DSD_Ref(manager, sub);
            __DSD_RecursiveDeref(manager, iter->decomposition);
            iter->decomposition = DSD_Regular(sub);


            if(DSD_IsComplement(sub))
            {
                complemented = 1;
                complemented_var = canonical_var(iter->decomposition);
            }


        }

        iter = iter->next;
    }

    result = create_DSD_node(manager, f);

    __DSD_Ref(manager, result);
    
    result->type_actualsize = node_reg->type_actualsize;
    SET_CAN((DSD_Regular(result)), (canonical_var(node_reg)));

    manager->total_actualsize += INPUT_SIZE(DSD_Regular(result));

    result->actual_list = conflict_residue->next;

#ifndef DISABLE_SBDD
    result->symbolic_kernel = node_reg->symbolic_kernel;
    Cudd_Ref(result->symbolic_kernel);

    if(complemented)
    {
        temp = result->symbolic_kernel;
        result->symbolic_kernel = symbolic_merger(manager->Ddmanager_analogue, temp, Cudd_Not(Cudd_bddIthVar(manager->Ddmanager_analogue, complemented_var)), Cudd_bddIthVar(manager->Ddmanager_analogue, complemented_var));
        Cudd_Ref(result->symbolic_kernel);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);
    }
#endif

    if(DSD_IsComplement(node))
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        result = DSD_Not(result);
    }

    /*erase top actual node that is not used in the new DSD actual list*/	
    FixHeapFree(actual_malloc_ptr, conflict_residue);

    return result;

}

ActualNode *cofactor_equivalence(DSDManager *manager, DSDNode *container, DSDNode *child)
{
    int found;

    ActualNode *iter1, *iter2, *result;
    DdNode *main_container, *main_child, *temp_container, *temp_child, *stupid1, *stupid2; 

#ifdef DISABLE_GC
    int *length;
    main_container = DSD_Regular(container)->bdd_analogue;
    main_child = DSD_Regular(child)->bdd_analogue;
    length = (int*) malloc(sizeof(int));
#else
    main_container = DSD_Regular(container)->symbolic_kernel;
    main_child = DSD_Regular(child)->symbolic_kernel;
#endif

    
    if(DSD_IsComplement(container))
    {
        main_container = Cudd_Not(main_container);
    }

    if(DSD_IsComplement(child))
    {
        main_child = Cudd_Not(main_child);
    }

    found = 0;

    for(iter1 = DSD_Regular(container)->actual_list, iter2 = DSD_Regular(child)->actual_list; iter1 && iter2; iter1 = iter1->next, iter2 = iter2->next)
    {

#ifdef DISABLE_GC
        temp_container = Cudd_LargestCube(manager->Ddmanager_analogue, iter1->decomposition->bdd_analogue, length);
        Cudd_Ref(temp_container);

        temp_child = Cudd_LargestCube(manager->Ddmanager_analogue, Cudd_Not(iter2->decomposition->bdd_analogue), length);
        Cudd_Ref(temp_child);
        
        stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, main_container, temp_container);
        Cudd_Ref(stupid1);
        stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, main_child, temp_child);
        Cudd_Ref(stupid2);
#else
        temp_container = Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(iter1->decomposition));
        temp_child = Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(iter2->decomposition));
        
        stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, main_container, temp_container);
        Cudd_Ref(stupid1);
        stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, main_child, Cudd_Not(temp_child));
        Cudd_Ref(stupid2);
#endif

        if(stupid1 == stupid2)
        {


           
            Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid1);
            Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid2);
#ifdef DISABLE_GC
            Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_container);
            Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_child);

            temp_container = Cudd_LargestCube(manager->Ddmanager_analogue, Cudd_Not(iter1->decomposition->bdd_analogue), length);
            Cudd_Ref(temp_container);
            temp_child = Cudd_LargestCube(manager->Ddmanager_analogue, iter2->decomposition->bdd_analogue, length);
            Cudd_Ref(temp_child);

            stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, main_container, temp_container);
            Cudd_Ref(stupid1);
            stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, main_child, temp_child);
            Cudd_Ref(stupid2);
#else
            stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, main_container, Cudd_Not(temp_container));
            Cudd_Ref(stupid1);
            stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, main_child, temp_child);
            Cudd_Ref(stupid2);
#endif



            if(stupid1 == stupid2)	
            {
                mark(DSD_Regular(iter1->decomposition));
                found = 1;
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid1);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid2);
#ifdef DISABLE_GC
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_container);
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_child);
#endif

               break; //break for DISABLE_GC
            }
        }

        Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid1);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid2);
#ifdef DISABLE_GC
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_container);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_child);
#endif

    }


    
#ifdef DISABLE_GC
    free(length);
#endif

    if(!found)
    {
        iter1 = DSD_Regular(container)->actual_list;

        while(iter1)
        {
            unmark(DSD_Regular(iter1->decomposition));
            iter1 = iter1->next;
        }

        return NULL;
    }


    return copy_actual_list(DSD_Regular(container)->actual_list);

}


/*returns signed ActualNode for negative cofactor equivalence, marks the sub node for substitution*/
ActualNode *cofactor_container_node_equivalence(DSDManager *manager, DSDNode *container, DSDNode *child)
{
    int can_var;
    DdNode *result_decomposition;
    DdNode *main_symbolic, *stupid;	

    DdNode *cof_result;
    DdNode *temp_var;

    int found;
    int negative;

   
    ActualNode *iter, *result, *iter2, *answer;

    can_var = canonical_var(child);

#ifdef DISABLE_GC
    int *length;
    result_decomposition = DSD_Regular(child)->bdd_analogue;
    length = (int *) malloc(sizeof(int));
    
    if(DSD_IsComplement(child))
    {
        result_decomposition = Cudd_Not(result_decomposition);
    }
    
    if(DSD_IsComplement(container))
    {
        main_symbolic = Cudd_Not(DSD_Regular(container)->bdd_analogue);
    }
    else
    {
        main_symbolic = DSD_Regular(container)->bdd_analogue;
    }
#else
    result_decomposition = Cudd_bddIthVar(manager->Ddmanager_analogue, can_var);
    
    if(DSD_IsComplement(child))
    {
        result_decomposition = Cudd_Not(result_decomposition);
    }

    if(DSD_IsComplement(container))
    {
        main_symbolic = Cudd_Not(DSD_Regular(container)->symbolic_kernel);
    }
    else
    {
        main_symbolic = DSD_Regular(container)->symbolic_kernel;
    }

#endif 
 
    found = 0;
    negative = 0;

    answer = NULL;

    iter = DSD_Regular(container)->actual_list;

    while(iter)
    {
        assert(!DSD_IsComplement(iter));

        if(iter->decomposition != child)
        {
                 
#ifdef DISABLE_GC        
            temp_var = Cudd_LargestCube(manager->Ddmanager_analogue, iter->decomposition->bdd_analogue, length);
            Cudd_Ref(temp_var);
#else
            temp_var = Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(iter->decomposition));
#endif

            /*assert(!DSD_IsComplement(iter->decomposition));*/
            stupid = Cudd_Cofactor(manager->Ddmanager_analogue, main_symbolic, temp_var);
            Cudd_Ref(stupid);


 
            if(stupid == result_decomposition)
            {

                Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid);
#ifdef DISABLE_GC
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_var);
#endif
                found = 1;
                answer = iter;
                break;
            }
            else
            {
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid);
#ifdef DISABLE_GC              
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_var);
                temp_var = Cudd_LargestCube(manager->Ddmanager_analogue, Cudd_Not(iter->decomposition->bdd_analogue), length);
                Cudd_Ref(temp_var);
                stupid = Cudd_Cofactor(manager->Ddmanager_analogue, main_symbolic, temp_var);
#else
                stupid = Cudd_Cofactor(manager->Ddmanager_analogue, main_symbolic, (Cudd_Not(temp_var)));
#endif           
                Cudd_Ref(stupid);

                if(stupid == result_decomposition)
                {

                    Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid);
#ifdef DISABLE_GC
                    Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_var);
#endif
                    found = 1;
                    negative = 1;
                    answer = iter;
                    break;
                }
            }

            Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid);
#ifdef DISABLE_GC
            Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_var);
#endif

        }



        iter = iter->next;
    }




#ifdef DISABLE_GC
    free(length);
#endif
    
    if(!found)
    {
        return NULL;
    }



    iter = DSD_Regular(container)->actual_list;

    if(iter == answer)
    {
        mark(DSD_Regular(iter->decomposition));
    }



    result = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    result->decomposition = iter->decomposition;
    result->next = NULL;

    iter2 = result;
    iter = iter->next;


    while(iter)
    {
        if(iter == answer)
        {
            mark(DSD_Regular(iter->decomposition));
        }

        iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
        iter2 = iter2->next;

        iter2->decomposition = iter->decomposition;
        iter2->next = NULL;

        iter = iter->next;
    }	



    if(negative)	{
        return (ActualNode *) DSD_Not(result);
    }

    return result;	

}


ActualNode *cofactor_cofactor_equivalence(DSDManager *manager, DSDNode *container, DSDNode *child, int *switch_them)
{

    int anti, *size, found;
    DdNode *Esymbolic, *Tsymbolic, *Ecof, *Tcof, *stupid1, *stupid2;
    ActualNode *Ttemp, *Etemp, *result;
    DSDNode *temp_node;

#ifdef DISABLE_GC
    int *length;
    length = (int*) malloc(sizeof(int));
#endif

    size = (int*) malloc(sizeof(int));
    
    anti = 0;
    found = 0;

#ifdef DISABLE_GC

    if(DSD_IsComplement(container))
    {
        Esymbolic = Cudd_Not(DSD_Regular(container)->bdd_analogue);
    }
    else
    {
        Esymbolic = container->bdd_analogue;
    }

    if(DSD_IsComplement(child))
    {
        Tsymbolic = Cudd_Not(DSD_Regular(child)->bdd_analogue);
    }
    else
    {
        Tsymbolic = child->bdd_analogue;
    }

#else

    if(DSD_IsComplement(container))
    {
        Esymbolic = Cudd_Not(DSD_Regular(container)->symbolic_kernel);
    }
    else
    {
        Esymbolic = container->symbolic_kernel;
    }

    if(DSD_IsComplement(child))
    {
        Tsymbolic = Cudd_Not(DSD_Regular(child)->symbolic_kernel);
    }
    else
    {
        Tsymbolic = child->symbolic_kernel;
    }

#endif
    
    Etemp = list_residue(manager->Ddmanager_analogue, DSD_Regular(container)->actual_list, DSD_Regular(child)->actual_list, size);

    assert(*size == 1);


    Ttemp = list_residue(manager->Ddmanager_analogue, DSD_Regular(child)->actual_list, DSD_Regular(container)->actual_list, size);

    assert(*size == 1);

    free(size);

#ifdef DISABLE_GC
    Ecof = Cudd_LargestCube(manager->Ddmanager_analogue, Etemp->decomposition->bdd_analogue, length);
    Cudd_Ref(Ecof);

    Tcof = Cudd_LargestCube(manager->Ddmanager_analogue, Ttemp->decomposition->bdd_analogue, length);
    Cudd_Ref(Tcof);
#else
    Ecof = Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(Etemp->decomposition));
    assert(!DSD_IsComplement(Etemp->decomposition));

    Tcof = Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(Ttemp->decomposition));
    assert(!DSD_IsComplement(Etemp->decomposition));
#endif

    stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, Esymbolic, Ecof);
    Cudd_Ref(stupid1);

    stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, Tsymbolic, Tcof);
    Cudd_Ref(stupid2);



    
    if(stupid1 == stupid2)
    {
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid1);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid2);

#ifdef DISABLE_GC
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, Ecof);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, Tcof);


        Ecof = Cudd_LargestCube(manager->Ddmanager_analogue, Cudd_Not(Etemp->decomposition->bdd_analogue), length);
        Cudd_Ref(Ecof);

        Tcof = Cudd_LargestCube(manager->Ddmanager_analogue, Cudd_Not(Ttemp->decomposition->bdd_analogue), length);
        Cudd_Ref(Tcof);

        stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, Esymbolic, Ecof);
        Cudd_Ref(stupid1);
        stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, Tsymbolic, Tcof);
        Cudd_Ref(stupid2);
#else
        stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, Esymbolic, Cudd_Not(Ecof));
        Cudd_Ref(stupid1);
        stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, Tsymbolic, Cudd_Not(Tcof));
        Cudd_Ref(stupid2);
#endif
        
        if(stupid1 == stupid2)
        {
            found = 1;
        }
    }
    else
    {
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid1);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid2);

#ifdef DISABLE_GC
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, Ecof);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, Tcof);

        Ecof = Cudd_LargestCube(manager->Ddmanager_analogue, Etemp->decomposition->bdd_analogue, length);
        Cudd_Ref(Ecof);

        Tcof = Cudd_LargestCube(manager->Ddmanager_analogue, Cudd_Not(Ttemp->decomposition->bdd_analogue), length);
        Cudd_Ref(Tcof);

        stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, Esymbolic, Ecof);
        Cudd_Ref(stupid1);

        stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, Tsymbolic, Tcof);
        Cudd_Ref(stupid2);
#else  
        stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, Esymbolic, Ecof);
        Cudd_Ref(stupid1);

        stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, Tsymbolic, Cudd_Not(Tcof));
        Cudd_Ref(stupid2);
#endif

        
        if(stupid1 == stupid2)
        {
            Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid1);
            Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid2);
#ifdef DISABLE_GC            
            Cudd_RecursiveDeref(manager->Ddmanager_analogue, Ecof);
            Cudd_RecursiveDeref(manager->Ddmanager_analogue, Tcof);

            Ecof = Cudd_LargestCube(manager->Ddmanager_analogue, Cudd_Not(Etemp->decomposition->bdd_analogue), length);
            Cudd_Ref(Ecof);

            Tcof = Cudd_LargestCube(manager->Ddmanager_analogue, Ttemp->decomposition->bdd_analogue, length);
            Cudd_Ref(Tcof);

            stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, Esymbolic, Ecof);
            Cudd_Ref(stupid1);
            stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, Tsymbolic, Tcof);
            Cudd_Ref(stupid2);
#else
            stupid1 = Cudd_Cofactor(manager->Ddmanager_analogue, Esymbolic, Cudd_Not(Ecof));
            Cudd_Ref(stupid1);
            stupid2 = Cudd_Cofactor(manager->Ddmanager_analogue, Tsymbolic, Tcof);
            Cudd_Ref(stupid2);
#endif
            
            if(stupid1 == stupid2)
            {
                anti = 1;
                found = 1;
            }
        }


    }



    
    Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid1);
    Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid2);

#ifdef DISABLE_GC
    Cudd_RecursiveDeref(manager->Ddmanager_analogue, Ecof);
    Cudd_RecursiveDeref(manager->Ddmanager_analogue, Tcof);
    free(length);
#endif


    if(!found)
    {
        FixHeapFree(actual_malloc_ptr, Etemp);
        FixHeapFree(actual_malloc_ptr, Ttemp);
        return NULL;
    }

    if(Cudd_ReadPerm(manager->Ddmanager_analogue, canonical_var(Etemp->decomposition)) < Cudd_ReadPerm(manager->Ddmanager_analogue, canonical_var(Ttemp->decomposition)))
    {
        mark(DSD_Regular(Ttemp->decomposition));
        result = Etemp;
        result->next = copy_actual_list(DSD_Regular(child)->actual_list);
        FixHeapFree(actual_malloc_ptr, Ttemp);
        *switch_them = 1;
    }
    else
    {
        mark(DSD_Regular(Etemp->decomposition));
        result = Ttemp;
        result->next = copy_actual_list(DSD_Regular(container)->actual_list);
        FixHeapFree(actual_malloc_ptr, Etemp);
    }





    if(anti)
    {
        return (ActualNode *) DSD_Not(result);
    }

    return result;

}

int symbolic_finder_builder(DSDManager *manager, DSDNode *node1, DSDNode *node2)
{
    int found;
    ActualNode *iter;
    DdNode *result, *temp_result, *symbolic_smasher, top_func;

    found = 0;


    if(marked(DSD_Regular(node2)))
    {
        unmark_recursive(DSD_Regular(node2));

        return 1;

    }
    else if(GET_TYPE(DSD_Regular(node2)) == DSD_VAR)
    {
        return 0;
    }

    iter = DSD_Regular(node2)->actual_list;

    while(iter)
    {
        if(!symbolic_finder_builder(manager, node1, iter->decomposition))
        {
            return 0;
        }

        iter = iter->next;
    }

    return 1;


}


int symbolic_finder_builder_incomplete(DSDManager *manager, DSDNode *node1, DSDNode *node2)
{
    int found,found2;
    ActualNode *iter;
    DdNode *result, *temp_result, *top_func;
    int symbolic_smasher;


    found = 0;
    found2 = 0;

    if(!marked(DSD_Regular(node2)))
    {
        return 1;

    }
    else if(GET_TYPE(DSD_Regular(node2)) == DSD_VAR)
    {
        return 2;
    }

    iter = DSD_Regular(node2)->actual_list;


    while(iter)
    {

        symbolic_smasher = symbolic_finder_builder_incomplete(manager, node1, iter->decomposition);

        if(symbolic_smasher == 1)
        {
            found = 1;
        }
        else if(symbolic_smasher == 0)
        {
            return 0;
        }
        else
        {
            found2 = 1;
        }



        if(found && found2)
        {
            return 0;
        }

        iter = iter->next;
    }



    /*!!!if(DSD_IsComplement(node2))
      {
      return Cudd_Not(result);
      }
      else
      {*/

    if(found) return 1;
    return 2;
    /*}*/


}


/*global list persists through recursion and contains pushed eligible nodes, with path, downn tree*/
/*contain a global list of created bdds in the expansion*/
/*when cofactor found--create new bdd with substitution--decompose, cofactor bdd (no longer generalized) of top parent wrt old bdd, create new bdd, decompose*/
/*delete all global structures, deref intermediate bdds, deref new bdd with sub, ref new base along with copys of residue*/

DSDNode *cofactor_elem_BDN_BDD_PRIME_SUB_ELEM(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *node1, DSDNode *node2)
{
    DSDNode *node1_reg, *node2_reg, *result;
    ActualNode *actuals_node1, *actuals_node2, *actuals_iter, *conflict_residue;

    DdNode *symbolic_smasher, *node2_expansion, *node1_expansion, *temp, *temp_var, *stupid;	

    /*newBdds *symbolic_list;*/

    int found, negative, corrupt, done;
    int *length;

    found = 0;
    negative = 0;
    corrupt = 0;
    done = 0;

    node1_reg = DSD_Regular(node1);
    node2_reg = DSD_Regular(node2);

    actuals_node1 = node1_reg->actual_list;
    actuals_node2 = node2_reg->actual_list;




    /*	if(GET_TYPE(node2_reg) == DSD_PRIME)  || (INPUT_SIZE(node2_reg) > INPUT_SIZE(node1_reg)))
        {*/
    node2_expansion = node2_reg->bdd_analogue;
    if(DSD_IsComplement(node2))
    {
        node2_expansion = Cudd_Not(node2_expansion);
    }

    if(DSD_IsComplement(node1))
    {
        node1_expansion = Cudd_Not(node1_reg->bdd_analogue);
    }
    else
    {
        node1_expansion = node1_reg->bdd_analogue;
    }



    //Cudd_Ref(node2_expansion);

    /*	actuals_iter = actuals_node1;

        while(actuals_iter)
        {
        mark(actuals_iter->decomposition);
        actuals_iter = actuals_iter->next;
        }
     */

    mark_recursive(node1, 0);

    actuals_iter = actuals_node2;

    while(actuals_iter)
    {
        if(!symbolic_finder_builder(manager, node1, actuals_iter->decomposition))
        {

            unmark_recursive(node1);

            return NULL;
        }

        /*all prime actuals should be positive*/
        /*temp = node2_expansion;
          node2_expansion = symbolic_merger(manager->Ddmanager_analogue, temp, symbolic_smasher, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actuals_iter->decomposition)));
          Cudd_Ref(node2_expansion);
          Cudd_RecursiveDeref(manager->Ddmanager_analogue, symbolic_smasher);
          Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);			*/
        actuals_iter = actuals_iter->next;
    }

    actuals_iter = actuals_node1;


    while(actuals_iter)
    {
        if(!symbolic_finder_builder_incomplete(manager, node2, actuals_iter->decomposition))
        {
            unmark_recursive(node1);

            return NULL;
        }
        /*	else if(symbolic_smasher)
                {
        //all prime actuals should be positive
        temp = node1_expansion;
        node1_expansion = symbolic_merger(manager->Ddmanager_analogue, temp, symbolic_smasher, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actuals_iter->decomposition)));
        Cudd_Ref(node1_expansion);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, symbolic_smasher);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp);			
        }*/
        actuals_iter = actuals_iter->next;
    }



    actuals_iter = actuals_node1;

    length = (int*) malloc(sizeof(int));

    while(actuals_iter)
    {
        done = 0;

        if(marked(DSD_Regular(actuals_iter->decomposition)) && !found)
        {

            temp_var = Cudd_LargestCube(manager->Ddmanager_analogue, actuals_iter->decomposition->bdd_analogue, length);
            Cudd_Ref(temp_var);

            /*temp_var = Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(actuals_iter->decomposition));*/

            assert(!DSD_IsComplement(actuals_iter->decomposition));

            stupid = Cudd_Cofactor(manager->Ddmanager_analogue, node1_expansion, temp_var);
            Cudd_Ref(stupid);

            if(stupid == node2_expansion)
            {
                //Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid);

                found = 1;
                done = 1;
            }
            else
            {
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_var);
                temp_var = Cudd_LargestCube(manager->Ddmanager_analogue, Cudd_Not(actuals_iter->decomposition->bdd_analogue), length);
                Cudd_Ref(temp_var);

                Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid);

                stupid = Cudd_Cofactor(manager->Ddmanager_analogue, node1_expansion, temp_var);
                Cudd_Ref(stupid);

                if(stupid == node2_expansion)
                {
                    //Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid);

                    found = 1;
                    negative = 1;
                    done = 1;
                }

            }

            Cudd_RecursiveDeref(manager->Ddmanager_analogue, stupid);
            Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_var);


        }



        unmark_recursive(actuals_iter->decomposition);

        if(done) mark(DSD_Regular(actuals_iter->decomposition));
        actuals_iter = actuals_iter->next;
    }

    free(length);





    unmark(DSD_Regular(node1));


    if(!found)
    {
        return NULL;
    }


    conflict_residue = copy_actual_list(actuals_node1);

    if(negative)
    {
        conflict_residue = (ActualNode *) DSD_Not(conflict_residue);
    }

    return BDN_BDD_PRIME_SUB(manager, f, top_func, conflict_residue, node1);



    return NULL;	
}


DdNode *check_symbolic2(DSDManager *manager, DSDNode *node)
{
    int found;
    ActualNode *iter;
    DdNode *result, *temp_result, *symbolic_smasher, top_func;

    found = 0;


    if(GET_TYPE(DSD_Regular(node)) == DSD_VAR)
    {
        assert(DSD_Regular(node)->bdd_analogue == DSD_Regular(node)->symbolic_kernel);

        return DSD_Regular(node)->bdd_analogue;
    }

    iter = DSD_Regular(node)->actual_list;

    result = DSD_Regular(node)->symbolic_kernel;

    while(iter)
    {
        symbolic_smasher = check_symbolic2(manager, iter->decomposition);

        if(GET_TYPE(DSD_Regular(node)) != DSD_OR)
        {
            assert(!DSD_IsComplement(iter->decomposition));
        }

        if(GET_TYPE(DSD_Regular(node)) == DSD_OR)
        {
            assert((GET_TYPE(DSD_Regular(iter->decomposition)) != DSD_OR) || DSD_IsComplement(iter->decomposition));
        }

        if(GET_TYPE(DSD_Regular(node)) == DSD_XOR)
        {
            assert((GET_TYPE(iter->decomposition) != DSD_XOR));
        }


        /*all prime actuals should be positive*/
        temp_result = result;


        result = symbolic_merger(manager->Ddmanager_analogue, temp_result, symbolic_smasher, Cudd_bddIthVar(manager->Ddmanager_analogue, canonical_var(iter->decomposition)));

        Cudd_Ref(result);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, symbolic_smasher);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, temp_result);	

        iter = iter->next;
    }

    assert(result == DSD_Regular(node)->bdd_analogue);

    return result;

}

void check_symbolic(DSDNode *node)
{
    ActualNode *iter;
    int not_var;

    not_var = 0;

    iter = DSD_Regular(node)->actual_list;

    while(iter)
    {
        if(GET_TYPE(DSD_Regular(iter->decomposition)) != DSD_VAR)
        {
            not_var = 1;
        }


    }


    if(!not_var)
    {
        assert((DSD_Regular(node)->symbolic_kernel) == (DSD_Regular(node)->bdd_analogue));
    } 

}






DSDNode* Prime_Decomp(DSDManager* manager, DdNode* f, DdNode* top_func, DSDNode* T, DSDNode* E)
{
    DSDNode *TReg;
    DSDNode *EReg;
    DdNode *one;
    DdNode *zero;

    DSDNode *Tresidue_dec;
    DSDNode *Eresidue_dec;
    DSDNode *result;
    
    DSDNode *mux_dec;



    ActualNode* commons;
    ActualNode* Tresidue;
    ActualNode* Eresidue;

    ActualNode* temp;

    int *size, *switch_them;
    int size2;



    TReg = DSD_Regular(T);
    EReg = DSD_Regular(E);

    ActualNode* conflicting_term_residue;

    one = Cudd_ReadOne(manager->Ddmanager_analogue);
    zero = Cudd_Not(one);

    if(GET_TYPE(TReg) != DSD_PRIME && GET_TYPE(EReg) != DSD_PRIME)
    {
#ifdef DEBUG_PRINT
        printf("-------Section 3.1----------\n");
#endif

        return NULL;
    }


    if(GET_TYPE(EReg) == DSD_PRIME && node_exists(EReg->actual_list, TReg))
    {
        /*seperate calculation for negative prime*/
        if((conflicting_term_residue = cofactor_container_node_equivalence(manager, E, T)))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 3.2.1----------\n");
#endif

            result =  BDN_BDD_PRIME_SUB(manager, f, Cudd_Not(top_func), conflicting_term_residue, E); 
            if(result)
                DSD_Regular(result)->topvar_refsize = ((DSD_Regular(result)->topvar_refsize) & 0xffff0000) | ((DSD_Regular(result)->topvar_refsize & 0x0000ffff) - 1);
            return result;

        }
        else
        {
#ifdef DEBUG_PRINT
            printf("-------Section 3.2.2----------\n");
#endif

            return NULL;
        }
    }


    if(GET_TYPE(TReg) == DSD_PRIME && node_exists(TReg->actual_list, EReg))
    {
        /*seperate calculation for negative prime*/
        if((conflicting_term_residue = cofactor_container_node_equivalence(manager, T, E)))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 3.3.1----------\n");
#endif

            result = BDN_BDD_PRIME_SUB(manager, f, top_func, conflicting_term_residue, T); 
            if(result)
                DSD_Regular(result)->topvar_refsize = ((DSD_Regular(result)->topvar_refsize) & 0xffff0000) | ((DSD_Regular(result)->topvar_refsize & 0x0000ffff) - 1);
            return result;

        
        }
        else
        {
#ifdef DEBUG_PRINT
            printf("-------Section 3.3.2----------\n");
#endif

            return NULL;
        }
    }


    /*if T is greater than E, cofactor expanded E to get T*/
    if(GET_TYPE(EReg) == DSD_PRIME && GET_TYPE(TReg) != DSD_VAR && (GET_TYPE(TReg) != DSD_PRIME || INPUT_SIZE(EReg) > INPUT_SIZE(TReg)))
    {
#ifdef DEBUG_PRINT
        printf("-------Section 3.4----------\n");
#endif

        result = cofactor_elem_BDN_BDD_PRIME_SUB_ELEM(manager, f, Cudd_Not(top_func), E, T);	
        if(result)
                DSD_Regular(result)->topvar_refsize = ((DSD_Regular(result)->topvar_refsize) & 0xffff0000) | ((DSD_Regular(result)->topvar_refsize & 0x0000ffff) - 1);
            return result;

        
    }

    if(GET_TYPE(TReg) == DSD_PRIME && GET_TYPE(EReg) != DSD_VAR && (GET_TYPE(EReg) != DSD_PRIME || INPUT_SIZE(TReg) > INPUT_SIZE(EReg)))
    {
#ifdef DEBUG_PRINT
        printf("-------Section 3.5----------\n");
#endif

        result = cofactor_elem_BDN_BDD_PRIME_SUB_ELEM(manager, f, top_func, T, E); 
        if(result)
                DSD_Regular(result)->topvar_refsize = ((DSD_Regular(result)->topvar_refsize) & 0xffff0000) | ((DSD_Regular(result)->topvar_refsize & 0x0000ffff) - 1);
            return result;

        
    }




    if(!((GET_TYPE(TReg) == DSD_PRIME) && (GET_TYPE(EReg) == DSD_PRIME)))
    {
#ifdef DEBUG_PRINT
        printf("-------Section 3.6----------\n");
#endif

        return NULL;
    }

    size = (int*) malloc(sizeof(int));

    list_intersection_special(manager->Ddmanager_analogue, TReg->actual_list, EReg->actual_list, size);

    
    if(*size == INPUT_SIZE(TReg))
    {
        free(size);

        if((conflicting_term_residue = cofactor_equivalence(manager, T, E)))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 3.7.1----------\n");
#endif
            
            /*positive if T cofactor is positive--else negative*/
            result = BDN_BDD_PRIME_XOR_SUB(manager, f, top_func, conflicting_term_residue, T);
            if(result)
                DSD_Regular(result)->topvar_refsize = ((DSD_Regular(result)->topvar_refsize) & 0xffff0000) | ((DSD_Regular(result)->topvar_refsize & 0x0000ffff) - 1);
            return result;
        
        }
        else
        {
#ifdef DEBUG_PRINT
            printf("-------Section 3.7.2----------\n");
#endif

            return NULL;
        }
    }


    /*finds supports and compares to see if different--could implement intersection*/
    if(INPUT_SIZE(TReg) == 1 && INPUT_SIZE(EReg) == 1) //&& support_compare(manager, E, T) == -2)
    {	
#ifdef DEBUG_PRINT
        printf("-------Section 3.8----------\n");
#endif

        free(size);
        return BDN_MUX_VAR_DEC_DEC(manager, f, top_func, E, T);
    }

    *size = (INPUT_SIZE(TReg) - *size);

    if(*size > 1)
    {
#ifdef DEBUG_PRINT
        printf("-------Section 3.9----------\n");
#endif

        free(size);
        return NULL;
    }

    switch_them = (int*) malloc(sizeof(int));
    *switch_them = 0;

    if((conflicting_term_residue = cofactor_cofactor_equivalence(manager, E, T, switch_them)))
    {
#ifdef DEBUG_PRINT
        printf("-------Section 3.10/3.11----------\n");
#endif



        if(*switch_them)
        {
            free(switch_them);
            result = BDN_BDD_PRIME_MUX_SUB(manager, f, Cudd_Not(top_func), conflicting_term_residue, T);
            
            if(result)
                DSD_Regular(result)->topvar_refsize = ((DSD_Regular(result)->topvar_refsize) & 0xffff0000) | ((DSD_Regular(result)->topvar_refsize & 0x0000ffff) - 1);
            return result;

        
        }
        else
        {
            free(switch_them);
            result = BDN_BDD_PRIME_MUX_SUB(manager, f, top_func, conflicting_term_residue, E);
            if(result)
                DSD_Regular(result)->topvar_refsize = ((DSD_Regular(result)->topvar_refsize) & 0xffff0000) | ((DSD_Regular(result)->topvar_refsize & 0x0000ffff) - 1);
            return result;

        }
    }

    free(switch_them);
    return NULL;


}

