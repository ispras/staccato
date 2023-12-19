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


#include "DSDOrDecompose.h"


DSDNode *BDN_OR_VAR_EXP(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter1;
    ActualNode *iter2;

    int count;

    count = 1;

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);

    result = create_DSD_node(manager, f);
    SET_CAN((DSD_Regular(result)), (canonical_var(base)));


    SET_TYPE(result, OR);


    top_node = create_var(manager, top_func);

    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    result->actual_list->decomposition = top_node;
    result->actual_list->next = NULL;

    iter2 = result->actual_list;
    iter1 = DSD_Regular(base)->actual_list;

    while(iter1 != NULL)
    {
        iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
        iter2 = iter2->next;

        iter2->decomposition = iter1->decomposition;
        iter2->next = NULL;

        __DSD_Ref(manager, iter2->decomposition);
        iter1 = iter1->next;
        count++;
    }

#ifndef DISABLE_SBDD            
    result->symbolic_kernel = symbolic_or(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif


    SET_SIZE(result, count);
    manager->total_actualsize += count;
    return result;

}









DSDNode *BDN_NOR_VAR_EXP(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter1;
    ActualNode *iter2;

    int count;

    count = 1;

    result = create_DSD_node(manager, f);
    SET_CAN((DSD_Regular(result)), (canonical_var(base)));

    SET_TYPE(result, OR);

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);

    top_node = create_var(manager, top_func);

    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    result->actual_list->decomposition = top_node;
    result->actual_list->next = NULL;

    iter2 = result->actual_list;
    iter1 = DSD_Regular(base)->actual_list;

    while(iter1 != NULL)
    {
        iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
        iter2 = iter2->next;

        iter2->decomposition = iter1->decomposition;
        iter2->next = NULL;

        __DSD_Ref(manager, iter2->decomposition);
        iter1 = iter1->next;
        count++;
    }

#ifndef DISABLE_SBDD            
    result->symbolic_kernel = symbolic_or(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);	
#endif

    SET_SIZE(result, count);
    manager->total_actualsize += count;

    result->bdd_analogue = Cudd_Not(result->bdd_analogue);

    return DSD_Not(result);
}


DSDNode *BDN_NOR_VAR_DEC(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter2;

    result = create_DSD_node(manager, f);
    SET_CAN((DSD_Regular(result)), (canonical_var(base)));

    SET_TYPE(result, OR);

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);

    top_node = create_var(manager, top_func);

    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    result->actual_list->decomposition = top_node;

    result->actual_list->next = NULL;

    iter2 = result->actual_list;

    iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
    iter2 = iter2->next;

    iter2->decomposition = base;

    iter2->next = NULL;

    __DSD_Ref(manager, iter2->decomposition);


#ifndef DISABLE_SBDD            
    result->symbolic_kernel = symbolic_or(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, 2);
    manager->total_actualsize += 2;


    result->bdd_analogue = Cudd_Not(result->bdd_analogue);

    return DSD_Not(result);
}


DSDNode *BDN_OR_VAR_DEC(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter2;


    result = create_DSD_node(manager, f);
    SET_CAN((DSD_Regular(result)), (canonical_var(base)));

    SET_TYPE(result, OR);

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);

    top_node = create_var(manager, top_func);

    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    result->actual_list->decomposition = top_node;

    result->actual_list->next = NULL;

    iter2 = result->actual_list;

    iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
    iter2 = iter2->next;

    iter2->decomposition = base;

    iter2->next = NULL;

    __DSD_Ref(manager, iter2->decomposition);


#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_or(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, 2);
    manager->total_actualsize += 2;

    return result;
}

DSDNode *BDN_OR_DEC_ACTUALS(DSDManager *manager, DdNode *f, DSDNode *node, ActualNode *actuals)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter1, *temp, *previous;
    ActualNode *iter2;

    int count, found;

    count = 1;
    found = 0;

    result = create_DSD_node(manager, f);

    SET_TYPE(result, OR);

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);

    temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
    temp->decomposition = node;


    __DSD_Ref(manager, node);

    result->actual_list = actuals;
    iter1 = actuals;
    previous = NULL;


    while(iter1 != NULL)
    {
        __DSD_Ref(manager, iter1->decomposition);


        if(!found && Cudd_ReadPerm(manager->Ddmanager_analogue, canonical_var(temp->decomposition)) < Cudd_ReadPerm(manager->Ddmanager_analogue, canonical_var(iter1->decomposition)))
        {
            found = 1;

            if(!previous)
            {
                result->actual_list = temp;
                temp->next = actuals;
            }
            else
            {
                previous->next = temp;
                temp->next = iter1;
            }


        }

        previous = iter1;
        iter1 = iter1->next;
        count++;
    }

    if(!found)
    {
        previous->next = temp;
        temp->next = NULL;
    }

#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_or(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    set_canonical_var(result);


    SET_SIZE(result, count);
    manager->total_actualsize += count;

    return result;

}


DSDNode *BDN_NOR_DEC_ACTUALS(DSDManager *manager, DdNode *f, DSDNode *node, ActualNode *actuals)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter1, *temp, *previous;
    ActualNode *iter2;

    int count, found;

    count = 1;
    found = 0;

    result = create_DSD_node(manager, f);

    SET_TYPE(result, OR);

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);

    __DSD_Ref(manager, node);

    temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
    temp->decomposition = node;


    result->actual_list = actuals;
    iter1 = actuals;
    previous = NULL;

    while(iter1 != NULL)
    {
        __DSD_Ref(manager, iter1->decomposition);


        if(!found && Cudd_ReadPerm(manager->Ddmanager_analogue, canonical_var(temp->decomposition)) < Cudd_ReadPerm(manager->Ddmanager_analogue, canonical_var(iter1->decomposition)))
        {
            found = 1;

            if(!previous)
            {
                result->actual_list = temp;
                temp->next = actuals;
            }
            else
            {
                previous->next = temp;
                temp->next = iter1;
            }


        }

        previous = iter1;
        iter1 = iter1->next;
        count++;
    }

    if(!found)
    {
        previous->next = temp;
        temp->next = NULL;
    }


#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_or(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, count);
    manager->total_actualsize += count;

    result->bdd_analogue = Cudd_Not(result->bdd_analogue);
    set_canonical_var(result);

    return DSD_Not(result);
}

DSDNode *BDN_OR_DEC_DEC(DSDManager *manager, DdNode *f, DSDNode *node1, DSDNode *node2)
{
    DSDNode *result;


    ActualNode *iter2;

    DSDNode *temp;
    int var1, var2, temp_var;

    result = create_DSD_node(manager, f);

    SET_TYPE(result, OR);



    var1 = canonical_var(node1);
    var2 = canonical_var(node2);

    if(Cudd_ReadPerm(manager->Ddmanager_analogue, var1) > Cudd_ReadPerm(manager->Ddmanager_analogue, var2))
    {
        temp = node1;
        node1 = node2;
        node2 = temp;

        temp_var = var1;
        var1 = var2;
        var2 = temp_var;		
    }

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);
    assert(var1 != var2);
    assert(var1 < 10000);
    assert(var2 < 10000);

    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    result->actual_list->decomposition = node1;

    result->actual_list->next = NULL;

    iter2 = result->actual_list;	
    __DSD_Ref(manager, iter2->decomposition);

    iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
    iter2 = iter2->next;

    iter2->decomposition = node2;

    iter2->next = NULL;

    __DSD_Ref(manager, iter2->decomposition);


#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_or(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, 2);
    manager->total_actualsize += 2;

    /*manager->Ddmanager_analogue, result->bdd_analogue;*/
    set_canonical_var(result);

    return result;
}

DSDNode *BDN_NOR_DEC_DEC(DSDManager *manager, DdNode *f, DSDNode *node1, DSDNode *node2)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter2;

    DSDNode *temp;
    int var1, var2, temp_var;

    result = create_DSD_node(manager, f);



    SET_TYPE(result, OR);



    var1 = canonical_var(node1);
    var2 = canonical_var(node2);

    if(Cudd_ReadPerm(manager->Ddmanager_analogue, var1) > Cudd_ReadPerm(manager->Ddmanager_analogue, var2))
    {
        temp = node1;
        node1 = node2;
        node2 = temp;

        temp_var = var1;
        var1 = var2;
        var2 = temp_var;		
    }


    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);
    assert(var1 != var2);
    assert(var1 < 10000);
    assert(var2 < 10000);

    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    result->actual_list->decomposition = node1;

    result->actual_list->next = NULL;

    iter2 = result->actual_list;	
    __DSD_Ref(manager, iter2->decomposition);

    iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
    iter2 = iter2->next;

    iter2->decomposition = node2;

    iter2->next = NULL;

    __DSD_Ref(manager, iter2->decomposition);


#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_or(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, 2);
    manager->total_actualsize += 2;

    result->bdd_analogue = Cudd_Not(result->bdd_analogue);

    set_canonical_var(result);
    return DSD_Not(result);
}



DSDNode *BDN_OR_ACTUALS(DSDManager *manager, DdNode *f, ActualNode *actuals)
{
    DSDNode *top_node;
    DSDNode *result, *last;

    ActualNode *iter1;
    ActualNode *iter2;

    int count;

    count = 0;

    result = create_DSD_node(manager, f);

    SET_TYPE(result, OR);


    result->actual_list = actuals;

    iter1 = actuals;

    while(iter1 != NULL)
    {
        __DSD_Ref(manager, iter1->decomposition);
        last = iter1->decomposition;
        iter1 = iter1->next;
        count++;
    }

    SET_CAN((DSD_Regular(result)), (canonical_var(last)));

#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_or(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, count);
    manager->total_actualsize += count;

    return result;


}

DSDNode *BDN_NOR_ACTUALS(DSDManager *manager, DdNode *f, ActualNode *actuals)
{
    DSDNode *top_node;
    DSDNode *result, *last;

    ActualNode *iter1;
    ActualNode *iter2;

    int count;

    count = 0;

    result = create_DSD_node(manager, f);

    SET_TYPE(result, OR);


    result->actual_list = actuals;

    iter1 = actuals;

    while(iter1 != NULL)
    {
        __DSD_Ref(manager, iter1->decomposition);
        last = iter1->decomposition;
        iter1 = iter1->next;
        count++;
    }

#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_or(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, count);
    manager->total_actualsize += count;

    SET_CAN((DSD_Regular(result)), (canonical_var(last)));
    result->bdd_analogue = Cudd_Not(result->bdd_analogue);

    return DSD_Not(result);

}




DSDNode *BDN_BDD_OR_RESIDUE(DSDManager *manager, ActualNode *residue)
{
    DdNode *f, *tmp;
    ActualNode *iter;

    DSDNode *result;


    if(residue == NULL)
    {
        return NULL;
    }
    else if(residue->next == NULL)
    {
        __DSD_Ref(manager, residue->decomposition);
        return residue->decomposition;
    }



    if(DSD_IsComplement(residue->decomposition))
    {
        f = Cudd_Not(DSD_Regular(residue->decomposition)->bdd_analogue);
    }
    else
    {
        f = residue->decomposition->bdd_analogue;
    }

    iter = residue->next;

    Cudd_Ref(f);


    while(iter)
    {
        assert(iter->decomposition != NULL);

        if(DSD_IsComplement(iter->decomposition))
        {
            tmp = Cudd_bddOr(manager->Ddmanager_analogue, Cudd_Not(DSD_Regular(iter->decomposition)->bdd_analogue),f);
        }
        else
        {
            tmp = Cudd_bddOr(manager->Ddmanager_analogue, iter->decomposition->bdd_analogue,f);
        }

        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);
        f = tmp;

        iter = iter->next;
    }


    if((result = find_DSD_node(manager, f)) == NULL)
    {
        result = BDN_OR_ACTUALS(manager, f, residue);
    }
   
    __DSD_Ref(manager, result); 
    
    Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);


    return result;

}

DSDNode *BDN_BDD_NOR_RESIDUE(DSDManager *manager, ActualNode *residue)
{
    DdNode *f, *tmp;
    ActualNode *iter;

    DSDNode *result;


    if(residue == NULL)
    {
        return NULL;
    }
    else if(residue->next == NULL)
    {
        if(GET_TYPE(DSD_Regular(residue->decomposition)) != OR || !DSD_IsComplement(residue->decomposition))
        {
            __DSD_Ref(manager, residue->decomposition);
            
            return DSD_Not(residue->decomposition);
        }
        else
        {
            return NULL;
        }
    }



    if(DSD_IsComplement(residue->decomposition))
    {
        f = Cudd_Not(DSD_Regular(residue->decomposition)->bdd_analogue);
    }
    else 
    {
        f = residue->decomposition->bdd_analogue;
    }

    iter = residue->next;

    Cudd_Ref(f);


    while(iter)
    {
        assert(iter->decomposition != NULL);

        if(DSD_IsComplement(iter->decomposition))
        {
            tmp = Cudd_bddOr(manager->Ddmanager_analogue, Cudd_Not(DSD_Regular(iter->decomposition)->bdd_analogue),f);
        }
        else
        {
            tmp = Cudd_bddOr(manager->Ddmanager_analogue, iter->decomposition->bdd_analogue,f);
        }

        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);
        f = tmp;

        iter = iter->next;
    }

    if(residue->next) f = Cudd_Not(f);

    if((result = find_DSD_node(manager, f)) == NULL)
    {
        result = BDN_NOR_ACTUALS(manager, f, residue);

    }
    
    __DSD_Ref(manager, result);


    Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);

    return result;
}



DSDNode *BDN_BDD_NOR_VAR_ACTUALS(DSDManager *manager, DdNode *top_func, ActualNode *residue)
{
    DdNode *f, *tmp;
    DSDNode *top_node;
    ActualNode *iter, *result_actuals;

    DSDNode *result;


    if(residue == NULL)
    {
        return NULL;
    }
    else if(residue->next == NULL)
    {
        if(GET_TYPE(DSD_Regular(residue->decomposition)) != OR  || !DSD_IsComplement(residue->decomposition))
        {
            __DSD_Ref(manager, residue->decomposition);            
            return DSD_Not(residue->decomposition);
        }
        else
        {
            return NULL;
        }
    }

    f = top_func;

    top_node = create_var(manager, top_func);

    result_actuals = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    result_actuals->decomposition = top_node;
    result_actuals->next = residue;



    iter = residue;

    Cudd_Ref(f);


    while(iter)
    {
        assert(iter->decomposition != NULL);

        if(DSD_IsComplement(iter->decomposition))
        {
            tmp = Cudd_bddOr(manager->Ddmanager_analogue, Cudd_Not(DSD_Regular(iter->decomposition)->bdd_analogue),f);
        }
        else
        {
            tmp = Cudd_bddOr(manager->Ddmanager_analogue, iter->decomposition->bdd_analogue,f);
        }

        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);
        f = tmp;

        iter = iter->next;
    }

    if((result = find_DSD_node(manager, Cudd_Not(f))) == NULL)
    {
        result = BDN_NOR_ACTUALS(manager, Cudd_Not(f), result_actuals);

    }
    else
    {
        /*delete allocated memory*/
    }

    __DSD_Ref(manager, result);

    Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);

    return result;
}











DSDNode *BDN_BDD_NOR_VAR_DEC(DSDManager *manager, DdNode *top_func, DSDNode *node)
{
    DSDNode *result, *top_node;

    DdNode *f;


    if(DSD_IsComplement(node))
    {
        f = Cudd_Not(DSD_Regular(node)->bdd_analogue);
    }
    else
    {
        f = node->bdd_analogue;
    }


    f = Cudd_bddOr(manager->Ddmanager_analogue, top_func, f);

    f = Cudd_Not(f); 

    Cudd_Ref(f);


    if((result = find_DSD_node(manager, f)) == NULL)
    {
        if(GET_TYPE(DSD_Regular(node)) != OR || DSD_IsComplement(node))
        {
            result = BDN_NOR_VAR_DEC(manager, f, top_func, node);
        }
        else
        {
            top_node = create_var(manager, top_func);


            result = BDN_NOR_DEC_ACTUALS(manager, f, top_node, copy_actual_list(node->actual_list));

        }

    }
    
    __DSD_Ref(manager, result);
    
    Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);

    return result;

}

DSDNode* OR_Decomp(DSDManager* manager, DdNode* f, DdNode* top_func, DSDNode* T, DSDNode* E)
{
    DSDNode *TReg;
    DSDNode *EReg;
    DdNode *one;
    DdNode *zero;

    DSDNode *Tresidue_dec;
    DSDNode *Eresidue_dec;
    DSDNode *result;


    DSDNode *mux_dec, *residue_dec;


    ActualNode *commons;
    ActualNode *Tresidue;
    ActualNode *Eresidue;

    ActualNode* temp;

    int *size;

    TReg = DSD_Regular(T);
    EReg = DSD_Regular(E);

    one = Cudd_ReadOne(manager->Ddmanager_analogue);
    zero = Cudd_Not(one);



    if(((TReg->bdd_analogue == zero) && DSD_IsComplement(T)) || ((TReg->bdd_analogue == one) && (!DSD_IsComplement(T))))
    {
        /*-nor is now just or as it should be*/
        if(GET_TYPE(EReg) == OR && (!DSD_IsComplement(E)))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.1.1----------\n");
#endif

            return BDN_OR_VAR_EXP(manager, f, top_func, E);
        }
        else
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.1.2----------\n");
#endif

            return BDN_OR_VAR_DEC(manager, f, top_func, E);
        }
    }




    if(((EReg->bdd_analogue == zero) && DSD_IsComplement(E)) || ((EReg->bdd_analogue == one) && (!DSD_IsComplement(E))))
    {
        /*-nor is now just or as it should be*/
        if((GET_TYPE(TReg) == OR) && (!DSD_IsComplement(T)))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.2.1----------\n");
#endif

            return BDN_OR_VAR_EXP(manager, f, Cudd_Not(top_func), T);
        }
        else
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.2.2----------\n");
#endif

            return BDN_OR_VAR_DEC(manager, f, Cudd_Not(top_func), T);
        }
    }



    if(((EReg->bdd_analogue == one) && DSD_IsComplement(E)) || ((EReg->bdd_analogue == zero) && (!DSD_IsComplement(E))))
    {
        /*-nor is now just or as it should be*/
        if((GET_TYPE(TReg) == OR) && (DSD_IsComplement(T)))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.3.1----------\n");
#endif

            return BDN_NOR_VAR_EXP(manager, f, Cudd_Not(top_func), T);
        }
        else
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.3.2----------\n");
#endif

            return BDN_NOR_VAR_DEC(manager, f, Cudd_Not(top_func), DSD_Not(T));
        }
    }	


    if(((TReg->bdd_analogue == one) && DSD_IsComplement(T)) || ((TReg->bdd_analogue ==zero) && (!DSD_IsComplement(T))))
    {
        /*-nor is now just or as it should be*/
        if((GET_TYPE(EReg) == OR) && (DSD_IsComplement(E)))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.1.1opposite----------\n");
#endif

            return BDN_NOR_VAR_EXP(manager, f, top_func, E);
        }
        else
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.1.2opposite----------\n");
#endif

            return BDN_NOR_VAR_DEC(manager, f, top_func, DSD_Not(E));
        }
    }

    if(GET_TYPE(TReg) != OR && GET_TYPE(EReg) != OR)
    {
#ifdef DEBUG_PRINT
        printf("-------Section 1.4----------\n");
#endif

        return NULL;
    }	


    size = (int*) malloc(sizeof(int));

    if((GET_TYPE(TReg) == OR && (!DSD_IsComplement(T))) && (GET_TYPE(EReg) == OR && (!DSD_IsComplement(E))))
    {
        commons = list_intersection(manager->Ddmanager_analogue, TReg->actual_list, EReg->actual_list, size); 

        if(commons != NULL)
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.5.1a----------\n");
#endif

            Tresidue = list_residue(manager->Ddmanager_analogue, TReg->actual_list, commons, size);
            Eresidue = list_residue(manager->Ddmanager_analogue, EReg->actual_list, commons, size);

            Tresidue_dec = BDN_BDD_OR_RESIDUE(manager, Tresidue); /*returns null if Tresidue is null--list should not be destroyed*/
            Eresidue_dec = BDN_BDD_OR_RESIDUE(manager, Eresidue);

            free(size);


            assert(Tresidue_dec || Eresidue_dec);

            protect(manager, commons);

            mux_dec = BDN_BDD_MUX_VAR_DEC_DEC(manager, top_func, Eresidue_dec, Tresidue_dec);

             __DSD_RecursiveDeref(manager, Tresidue_dec);
             __DSD_RecursiveDeref(manager, Eresidue_dec);
            
            result =  BDN_OR_DEC_ACTUALS(manager, f, mux_dec, commons); 
        
            unprotect(manager, DSD_Regular(result)->actual_list);            
            return result;

        }


        free(size);

#ifdef DEBUG_PRINT
        printf("-------Section 1.5.2a----------\n");
#endif

        return NULL;	

    }



    /*both are nor's*/
    if((GET_TYPE(TReg) == OR && (DSD_IsComplement(T))) && (GET_TYPE(EReg) == OR && (DSD_IsComplement(E))))
    {
        commons = list_intersection(manager->Ddmanager_analogue, TReg->actual_list, EReg->actual_list, size); 

        if(commons != NULL)
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.5.1b----------\n");
#endif

            Tresidue = list_residue(manager->Ddmanager_analogue, TReg->actual_list, commons, size);
            Eresidue = list_residue(manager->Ddmanager_analogue, EReg->actual_list, commons, size);

            Tresidue_dec = BDN_BDD_OR_RESIDUE(manager, Tresidue); /*returns null if Tresidue is null--list should not be destroyed*/
            Eresidue_dec = BDN_BDD_OR_RESIDUE(manager, Eresidue);

            free(size);

            assert(Tresidue_dec || Eresidue_dec);

            protect(manager, commons);

            mux_dec = BDN_BDD_MUX_VAR_DEC_DEC(manager, top_func, Eresidue_dec, Tresidue_dec);

             __DSD_RecursiveDeref(manager, Tresidue_dec);
             __DSD_RecursiveDeref(manager, Eresidue_dec);
            
            result = BDN_NOR_DEC_ACTUALS(manager, f, mux_dec, commons); 
            
            unprotect(manager, DSD_Regular(result)->actual_list);
        
            return result;
        }


        free(size);	

#ifdef DEBUG_PRINT
        printf("-------Section 1.5.2b----------\n");
#endif

        return NULL;	

    }


    if(GET_TYPE(EReg) == OR)
    {

        if(node_exists(EReg->actual_list, T) && !DSD_IsComplement(E))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.6.1a----------\n");
#endif

            temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
            temp->next = NULL;
            temp->decomposition = T;

            /*not really a residue--diff between c0exp and c1*/
            Eresidue = list_residue(manager->Ddmanager_analogue, EReg->actual_list, temp, size);
            Eresidue_dec = BDN_BDD_NOR_RESIDUE(manager, Eresidue);

            if(Eresidue_dec)
            {			
                residue_dec = BDN_BDD_NOR_VAR_DEC(manager, top_func, Eresidue_dec);		
            }
            else
            {
                residue_dec = BDN_BDD_NOR_VAR_ACTUALS(manager, top_func, copy_actual_list(DSD_Regular(Eresidue->decomposition)->actual_list));
            }
           

            __DSD_RecursiveDeref(manager, Eresidue_dec);
            

            FixHeapFree(actual_malloc_ptr, temp);
            free(size);

            result = BDN_OR_DEC_DEC(manager, f, T, residue_dec);
            __DSD_RecursiveDeref(manager, residue_dec);
            return result;


        }
        else if(node_exists(EReg->actual_list, DSD_Not(T)) && DSD_IsComplement(E))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.6.1b----------\n");
#endif

            temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
            temp->next = NULL;
            temp->decomposition = DSD_Not(T);


            /*not really a residue--diff between c0exp and c1*/
            Eresidue = list_residue(manager->Ddmanager_analogue, EReg->actual_list, temp, size);
            Eresidue_dec = BDN_BDD_NOR_RESIDUE(manager, Eresidue);

            if(Eresidue_dec)
            {			
                residue_dec = BDN_BDD_NOR_VAR_DEC(manager, top_func, Eresidue_dec);		
            }
            else
            {
                residue_dec = BDN_BDD_NOR_VAR_ACTUALS(manager, top_func, copy_actual_list(DSD_Regular(Eresidue->decomposition)->actual_list));
            }

            __DSD_RecursiveDeref(manager, Eresidue_dec);
            
            free(size);
            FixHeapFree(actual_malloc_ptr, temp);

            result = BDN_NOR_DEC_DEC(manager, f, DSD_Not(T), residue_dec);
            __DSD_RecursiveDeref(manager, residue_dec);
            return result;

        }

#ifdef DEBUG_PRINT
        printf("-------Section 1.6.2----------\n");
#endif

    }


    if(GET_TYPE(TReg) == OR)
    {


        if(node_exists(TReg->actual_list, E) && !DSD_IsComplement(T))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.7.1a----------\n");
#endif

            temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
            temp->next = NULL;
            temp->decomposition = E;


            /*not really a residue--diff between c0exp and c1*/
            Tresidue = list_residue(manager->Ddmanager_analogue, TReg->actual_list, temp, size);
            Tresidue_dec = BDN_BDD_NOR_RESIDUE(manager, Tresidue);

            if(Tresidue_dec)
            {			
                residue_dec = BDN_BDD_NOR_VAR_DEC(manager, Cudd_Not(top_func), Tresidue_dec);		
            }
            else
            {
                residue_dec = BDN_BDD_NOR_VAR_ACTUALS(manager, Cudd_Not(top_func), copy_actual_list(DSD_Regular(Tresidue->decomposition)->actual_list));
            }		


            __DSD_RecursiveDeref(manager, Tresidue_dec);

            free(size);
            FixHeapFree(actual_malloc_ptr, temp);

            result = BDN_OR_DEC_DEC(manager, f, E, residue_dec);
            __DSD_RecursiveDeref(manager, residue_dec);        
            return result;


        }
        else if(node_exists(TReg->actual_list, DSD_Not(E)) && DSD_IsComplement(T))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 1.7.1b----------\n");
#endif			

            temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
            temp->next = NULL;
            temp->decomposition = DSD_Not(E);


            /*not really a residue--diff between c0exp and c1*/
            Tresidue = list_residue(manager->Ddmanager_analogue, TReg->actual_list, temp, size);
            Tresidue_dec = BDN_BDD_NOR_RESIDUE(manager, Tresidue);

            if(Tresidue_dec)
            {			
                residue_dec = BDN_BDD_NOR_VAR_DEC(manager, Cudd_Not(top_func), Tresidue_dec);		
            }
            else
            {
                residue_dec = BDN_BDD_NOR_VAR_ACTUALS(manager, Cudd_Not(top_func), copy_actual_list(DSD_Regular(Tresidue->decomposition)->actual_list));
            }

            __DSD_RecursiveDeref(manager, Tresidue_dec);
            
            free(size);
            FixHeapFree(actual_malloc_ptr, temp);

            result = BDN_NOR_DEC_DEC(manager, f, DSD_Not(E), residue_dec);
            __DSD_RecursiveDeref(manager, residue_dec);        
            return result;
        }

#ifdef DEBUG_PRINT
        printf("-------Section 1.7.2----------\n");
#endif

    }

    return NULL;

}

