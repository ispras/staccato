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


#include "DSDXorDecompose.h"


DSDNode *BDN_NXOR_VAR_EXP(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter1;
    ActualNode *iter2;

    int count;

    int negative;

    count = 1;

    negative = 1;

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);

    result = create_DSD_node(manager, f);
    SET_CAN((DSD_Regular(result)), (canonical_var(base)));

    SET_TYPE(result, DSD_XOR);


    top_node = create_var(manager, top_func);


    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    if(Cudd_IsComplement(top_node))
    {
        result->actual_list->decomposition = DSD_Not(top_node);
        negative++;
    }
    else
    {
        result->actual_list->decomposition = top_node;
    }
    result->actual_list->decomposition = top_node;
    result->actual_list->next = NULL;

    iter2 = result->actual_list;
    iter1 = DSD_Regular(base)->actual_list;

    while(iter1 != NULL)
    {
        iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
        iter2 = iter2->next;

        if(DSD_IsComplement(iter1->decomposition))
        {
            iter2->decomposition = DSD_Not(iter1->decomposition);
            negative++;
        }
        else
        {
            iter2->decomposition = iter1->decomposition;
        }

        iter2->next = NULL;

        __DSD_Ref(manager, iter2->decomposition);
        iter1 = iter1->next;
        count++;
    }

#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_xor(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, count);
    manager->total_actualsize += count;

    if((negative % 2) == 1)
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        result = DSD_Not(result);
    }

    return result;
}


DSDNode *BDN_XOR_VAR_EXP(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter1;
    ActualNode *iter2;

    int count;

    int negative;

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);

    count = 1;

    negative = 0;

    result = create_DSD_node(manager, f);

    SET_CAN((DSD_Regular(result)), (canonical_var(base)));


    SET_TYPE(result, DSD_XOR);


    top_node = create_var(manager, top_func);

    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    if(Cudd_IsComplement(top_node))
    {
        result->actual_list->decomposition = DSD_Not(top_node);
        negative++;
    }
    else
    {
        result->actual_list->decomposition = top_node;
    }

    result->actual_list->next = NULL;

    iter2 = result->actual_list;
    iter1 = DSD_Regular(base)->actual_list;

    while(iter1 != NULL)
    {
        iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
        iter2 = iter2->next;

        if(DSD_IsComplement(iter1->decomposition))
        {
            iter2->decomposition = DSD_Not(iter1->decomposition);
            negative++;
        }
        else
        {
            iter2->decomposition = iter1->decomposition;
        }


        iter2->next = NULL;

        __DSD_Ref(manager, iter2->decomposition);
        iter1 = iter1->next;
        count++;
    }

#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_xor(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, count);
    manager->total_actualsize += count;

    if((negative % 2) == 1)
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        result = DSD_Not(result);
    }

    return result;

}

DSDNode *BDN_NXOR_VAR_DEC(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter2;
    int negative = 1;

    result = create_DSD_node(manager, f);

    SET_CAN((DSD_Regular(result)), (canonical_var(base)));

    SET_TYPE(result, DSD_XOR);

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);

    top_node = create_var(manager, top_func);



    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    if(DSD_IsComplement(top_node))
    {
        result->actual_list->decomposition = DSD_Not(top_node);
        negative++;
    }
    else
    {
        result->actual_list->decomposition = top_node;
    }


    result->actual_list->next = NULL;

    iter2 = result->actual_list;

    iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
    iter2 = iter2->next;

    if(DSD_IsComplement(base))
    {
        iter2->decomposition = DSD_Not(base);
        negative++;
    }
    else
    {
        iter2->decomposition = base;
    }


    iter2->next = NULL;

    __DSD_Ref(manager, iter2->decomposition);


#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_xor(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, 2);
    manager->total_actualsize += 2;

    if(negative % 2)
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        return DSD_Not(result);
    }


    return result;
}


DSDNode *BDN_XOR_VAR_DEC(DSDManager *manager, DdNode *f, DdNode *top_func, DSDNode *base)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter2;
    int negative = 0;

    result = create_DSD_node(manager, f);

    SET_CAN((DSD_Regular(result)), (canonical_var(base)));


    SET_TYPE(result, DSD_XOR);



    top_node = create_var(manager, top_func);



    result->actual_list = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    if(DSD_IsComplement(top_node))
    {
        result->actual_list->decomposition = DSD_Not(top_node);
        negative++;
    }
    else
    {
        result->actual_list->decomposition = top_node;
    }


    result->actual_list->next = NULL;

    iter2 = result->actual_list;

    iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
    iter2 = iter2->next;

    if(DSD_IsComplement(base))
    {
        iter2->decomposition = DSD_Not(base);
        negative++;
    }
    else
    {
        iter2->decomposition = base;
    }


    iter2->next = NULL;

    __DSD_Ref(manager, iter2->decomposition);

#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_xor(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, 2);
    manager->total_actualsize += 2;

    if(negative % 2)
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        return DSD_Not(result);
    }


    return result;
}


DSDNode *BDN_XOR_DEC_ACTUALS(DSDManager *manager, DdNode *f, DSDNode *node, ActualNode *actuals)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter1, *temp, *previous;
    ActualNode *iter2;

    int count, found;

    int negative;

    count = 1;
    found = 0;

    negative = 0;

    result = create_DSD_node(manager, f);

    SET_TYPE(result, DSD_XOR);

    __DSD_Ref(manager, node);

    temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    if(DSD_IsComplement(node))
    {
        temp->decomposition = DSD_Not(node);
        negative++;
    }
    else
    {
        temp->decomposition = node;
    }	

    result->actual_list = actuals;
    iter1 = actuals;
    previous = NULL;


    while(iter1 != NULL)
    {

        if(DSD_IsComplement(iter1->decomposition))
        {
            iter1->decomposition = DSD_Not(iter1->decomposition);
            negative++;
        }


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
    result->symbolic_kernel = symbolic_xor(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, count);
    manager->total_actualsize += count;

    if((negative % 2) == 1)
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        result = DSD_Not(result);
    }
    set_canonical_var(result);

    return result;
}

DSDNode *BDN_NXOR_DEC_ACTUALS(DSDManager *manager, DdNode *f, DSDNode *node, ActualNode *actuals)
{
    DSDNode *result;

    ActualNode *iter1, *temp, *previous;
    ActualNode *iter2;

    int count, found;

    int negative;

    count = 1;
    found = 0;

    negative = 1;

    __DSD_Ref(manager, node);

    result = create_DSD_node(manager, f);

    SET_TYPE(result, DSD_XOR);

    assert(Cudd_Regular(f) != DD_ONE(manager->Ddmanager_analogue));
    assert(f != NULL);	


    temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);

    if(DSD_IsComplement(node))
    {
        temp->decomposition = DSD_Not(node);
        negative++;
    }
    else
    {
        temp->decomposition = node;
    }	

    result->actual_list = actuals;
    iter1 = actuals;
    previous = NULL;


    while(iter1 != NULL)
    {

        if(DSD_IsComplement(iter1->decomposition))
        {
            iter1->decomposition = DSD_Not(iter1->decomposition);
            negative++;
        }

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
    result->symbolic_kernel = symbolic_xor(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, count);
    manager->total_actualsize += count;

    if((negative % 2) == 1)
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        result = DSD_Not(result);
    }

    set_canonical_var(result);
    return result;
}



DSDNode *BDN_XOR_DEC_DEC(DSDManager *manager, DdNode *f, DSDNode *node1, DSDNode *node2)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter2;

    DSDNode *temp;
    int var1, var2, temp_var;

    int negative = 0;

    result = create_DSD_node(manager, f);



    SET_TYPE(result, DSD_XOR);



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

    if(DSD_IsComplement(node1))
    {
        result->actual_list->decomposition = DSD_Not(node1);
        negative++;
    }
    else
    {
        result->actual_list->decomposition = node1;
    }


    result->actual_list->next = NULL;

    iter2 = result->actual_list;

    __DSD_Ref(manager, iter2->decomposition);

    iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
    iter2 = iter2->next;

    if(DSD_IsComplement(node2))
    {
        iter2->decomposition = DSD_Not(node2);
        negative++;
    }
    else
    {
        iter2->decomposition = node2;
    }


    iter2->next = NULL;

    __DSD_Ref(manager, iter2->decomposition);


#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_xor(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, 2);
    manager->total_actualsize += 2;

    set_canonical_var(result);

    if(negative % 2)
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        return DSD_Not(result);
    }


    return result;
}


DSDNode *BDN_NXOR_DEC_DEC(DSDManager *manager, DdNode *f, DSDNode *node1, DSDNode *node2)
{
    DSDNode *top_node;
    DSDNode *result;

    ActualNode *iter2;

    DSDNode *temp;
    int var1, var2, temp_var;
    int negative = 1;

    result = create_DSD_node(manager, f);



    SET_TYPE(result, DSD_XOR);



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

    if(DSD_IsComplement(node1))
    {
        result->actual_list->decomposition = DSD_Not(node1);
        negative++;
    }
    else
    {
        result->actual_list->decomposition = node1;
    }


    result->actual_list->next = NULL;

    iter2 = result->actual_list;

    __DSD_Ref(manager, iter2->decomposition);

    iter2->next = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
    iter2 = iter2->next;

    if(DSD_IsComplement(node2))
    {
        iter2->decomposition = DSD_Not(node2);
        negative++;
    }
    else
    {
        iter2->decomposition = node2;
    }


    iter2->next = NULL;

    __DSD_Ref(manager, iter2->decomposition);


#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_xor(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_SIZE(result, 2);
    manager->total_actualsize += 2;
    set_canonical_var(result);

    if(negative % 2)
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        return DSD_Not(result);
    }


    return result;
}


DSDNode *BDN_XOR_ACTUALS(DSDManager *manager, DdNode *f, ActualNode *actuals)
{
    DSDNode *top_node;
    DSDNode *result, *last;

    ActualNode *iter1;
    ActualNode *iter2;

    int count;

    int negative;

    count = 0;

    negative = 0;

    result = create_DSD_node(manager, f);

    SET_TYPE(result, DSD_XOR);



    result->actual_list = actuals;

    iter1 = actuals;

    while(iter1 != NULL)
    {

        if(DSD_IsComplement(iter1->decomposition))
        {
            iter1->decomposition = DSD_Not(iter1->decomposition);
            negative++;
        }
        else
        {
            iter1->decomposition = iter1->decomposition;
        }

        last = iter1->decomposition;
        __DSD_Ref(manager, iter1->decomposition);
        iter1 = iter1->next;
        count++;
    }

#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_xor(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_CAN((DSD_Regular(result)), (canonical_var(last)));

    SET_SIZE(result, count);
    manager->total_actualsize += count;

    if((negative % 2) == 1)
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        result = DSD_Not(result);
    }

    return result;

}

DSDNode *BDN_NXOR_ACTUALS(DSDManager *manager, DdNode *f, ActualNode *actuals)
{
    DSDNode *top_node;
    DSDNode *result, *last;

    ActualNode *iter1;
    ActualNode *iter2;

    int count;

    int negative;

    count = 0;

    negative = 1;

    result = create_DSD_node(manager, f);

    SET_TYPE(result, DSD_XOR);


    result->actual_list= actuals;

    iter1 = actuals;

    while(iter1 != NULL)
    {

        if(DSD_IsComplement(iter1->decomposition))
        {
            iter1->decomposition = DSD_Not(iter1->decomposition);
            negative++;
        }
        else
        {
            iter1->decomposition = iter1->decomposition;
        }
        last = iter1->decomposition;

        __DSD_Ref(manager, iter1->decomposition);
        iter1 = iter1->next;
        count++;
    }

#ifndef DISABLE_SBDD
    result->symbolic_kernel = symbolic_xor(manager->Ddmanager_analogue, result);
    Cudd_Ref(result->symbolic_kernel);
#endif

    SET_CAN((DSD_Regular(result)), (canonical_var(last)));



    SET_SIZE(result, count);
    manager->total_actualsize += count;

    if((negative % 2) == 1)
    {
        result->bdd_analogue = Cudd_Not(result->bdd_analogue);
        result = DSD_Not(result);
    }

    return result;

}




DSDNode *BDN_BDD_XOR_RESIDUE(DSDManager *manager, ActualNode *residue)
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
            tmp = Cudd_bddXor(manager->Ddmanager_analogue, Cudd_Not(DSD_Regular(iter->decomposition)->bdd_analogue),f);
        }
        else
        {
            tmp = Cudd_bddXor(manager->Ddmanager_analogue, iter->decomposition->bdd_analogue,f);
        }

        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);
        f = tmp;

        iter = iter->next;
    }

    if((result = find_DSD_node(manager, f)) == NULL)
    {
        result = BDN_XOR_ACTUALS(manager, f, residue);
    }

    __DSD_Ref(manager, result);

    Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);


    return result;
}


DSDNode *BDN_BDD_NXOR_RESIDUE(DSDManager *manager, ActualNode *residue)
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
        return DSD_Not(residue->decomposition);
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
            tmp = Cudd_bddXor(manager->Ddmanager_analogue, Cudd_Not(DSD_Regular(iter->decomposition)->bdd_analogue),f);
        }
        else
        {
            tmp = Cudd_bddXor(manager->Ddmanager_analogue, iter->decomposition->bdd_analogue,f);
        }

        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);
        f = tmp;

        iter = iter->next;
    }

    if(residue->next) f = Cudd_Not(f);

    if((result = find_DSD_node(manager, f)) == NULL)
    {
        result = BDN_NXOR_ACTUALS(manager, f, residue);
    }
    __DSD_Ref(manager, result);

    Cudd_RecursiveDeref(manager->Ddmanager_analogue, f);


    return result;
}


DSDNode* XOR_Decomp(DSDManager* manager, DdNode* f, DdNode *top_func, DSDNode* T, DSDNode* E)
{
    DSDNode *TReg;
    DSDNode *EReg;
    DdNode *one;
    DdNode *zero;

    DSDNode *Tresidue_dec;
    DSDNode *Eresidue_dec;

    DSDNode *mux_dec;

    DSDNode *result;

    ActualNode* commons;
    ActualNode* Tresidue;
    ActualNode* Eresidue;

    ActualNode* temp;

    int *size;



    TReg = DSD_Regular(T);
    EReg = DSD_Regular(E);

    one = Cudd_ReadOne(manager->Ddmanager_analogue);
    zero = Cudd_Not(one);

    if(T == DSD_Not(E))
    {
        if(GET_TYPE(TReg) == DSD_XOR && (!DSD_IsComplement(T)))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 2.1.1----------\n");
#endif

            return BDN_NXOR_VAR_EXP(manager, f, top_func, T);
        }

        if(GET_TYPE(TReg) == DSD_XOR && (DSD_IsComplement(T)))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 2.1.2----------\n");
#endif

            return BDN_XOR_VAR_EXP(manager, f, top_func, T);
        }

#ifdef DEBUG_PRINT
        printf("-------Section 2.1.3----------\n");
#endif

        return BDN_NXOR_VAR_DEC(manager, f, top_func, T);

    }


    if(GET_TYPE(TReg) != DSD_XOR && GET_TYPE(EReg) != DSD_XOR)
    {
#ifdef DEBUG_PRINT
        printf("-------Section 2.2---------\n");
#endif

        return NULL;
    }



    size = (int*) malloc(sizeof(int));


    if((GET_TYPE(TReg) == DSD_XOR) && (GET_TYPE(EReg) == DSD_XOR))
    {
#ifdef DEBUG_PRINT
        printf("-------Section 2.3.1----------\n");
#endif



        commons = list_intersection(manager->Ddmanager_analogue, TReg->actual_list, EReg->actual_list, size);

        if(commons != NULL)
        {

            Tresidue = list_residue(manager->Ddmanager_analogue, TReg->actual_list, commons, size);
            Eresidue = list_residue(manager->Ddmanager_analogue, EReg->actual_list, commons, size);

            if DSD_IsComplement(T)
            {
                Tresidue_dec = BDN_BDD_NXOR_RESIDUE(manager, Tresidue);
                if(!Tresidue_dec)
                {
                    Tresidue_dec = manager->one;
                }				
            }
            else
                Tresidue_dec = BDN_BDD_XOR_RESIDUE(manager, Tresidue);

            if DSD_IsComplement(E)
            {
                Eresidue_dec = BDN_BDD_NXOR_RESIDUE(manager, Eresidue);
                if(!Eresidue_dec)
                {
                    Eresidue_dec = manager->one;
                }			
            }
            else
                Eresidue_dec = BDN_BDD_XOR_RESIDUE(manager, Eresidue);


            protect(manager, commons);

            mux_dec = BDN_BDD_MUX_VAR_DEC_DEC(manager, top_func, Eresidue_dec, Tresidue_dec);


            __DSD_RecursiveDeref(manager, Eresidue_dec);
            __DSD_RecursiveDeref(manager, Tresidue_dec);

            free(size);

            result =  BDN_XOR_DEC_ACTUALS(manager, f, mux_dec, commons);

            unprotect(manager, DSD_Regular(result)->actual_list);

            return result;
        }

        free(size);

        return NULL;

    }



    if((GET_TYPE(EReg)) == DSD_XOR && node_exists(EReg->actual_list, TReg) && DSD_IsComplement(T))
    {
#ifdef DEBUG_PRINT
        printf("-------Section 2 NEW 1----------\n");
#endif

        temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
        temp->next = NULL;
        temp->decomposition = TReg;

        Eresidue = list_residue(manager->Ddmanager_analogue, EReg->actual_list, temp, size);

        Eresidue_dec = BDN_BDD_XOR_RESIDUE(manager, Eresidue);

        if(!DSD_IsComplement(E))
        {
            Eresidue_dec = DSD_Not(Eresidue_dec);
        }			

        mux_dec = BDN_BDD_NOR_VAR_DEC(manager, top_func, DSD_Not(Eresidue_dec));


        __DSD_RecursiveDeref(manager, Eresidue_dec);

        free(size);
        FixHeapFree(actual_malloc_ptr, temp);

        result = BDN_XOR_DEC_DEC(manager, f, mux_dec, T);
        __DSD_RecursiveDeref(manager, mux_dec);        
        return result;


    }

    if((GET_TYPE(TReg)) == DSD_XOR && node_exists(TReg->actual_list, EReg) && DSD_IsComplement(E))
    {
#ifdef DEBUG_PRINT
        printf("-------Section 2 NEW 2----------\n");
#endif

        temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
        temp->next = NULL;
        temp->decomposition = EReg;


        Tresidue = list_residue(manager->Ddmanager_analogue, TReg->actual_list, temp, size);

        Tresidue_dec = BDN_BDD_XOR_RESIDUE(manager, Tresidue);

        if(!DSD_IsComplement(T))
        {
            Tresidue_dec = DSD_Not(Tresidue_dec);
        }			

        mux_dec = BDN_BDD_NOR_VAR_DEC(manager, Cudd_Not(top_func), DSD_Not(Tresidue_dec));


        __DSD_RecursiveDeref(manager, Tresidue_dec);


        free(size);
        FixHeapFree(actual_malloc_ptr, temp);

        result = BDN_XOR_DEC_DEC(manager, f, mux_dec, E);
        __DSD_RecursiveDeref(manager, mux_dec);        
        return result;

    }




    if((GET_TYPE(EReg)) == DSD_XOR)
    {
        if(node_exists(EReg->actual_list, T))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 2.4.1----------\n");
#endif

            temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
            temp->next = NULL;
            temp->decomposition = T;

            Eresidue = list_residue(manager->Ddmanager_analogue, EReg->actual_list, temp, size);

            Eresidue_dec = BDN_BDD_XOR_RESIDUE(manager, Eresidue);

            if(DSD_IsComplement(E))
            {
                Eresidue_dec = DSD_Not(Eresidue_dec);
            }			

            mux_dec = BDN_BDD_MUX_VAR_DEC_DEC(manager, top_func, Eresidue_dec, DSD_Not(manager->one));

            __DSD_RecursiveDeref(manager, Eresidue_dec);

            free(size);
            FixHeapFree(actual_malloc_ptr, temp);

            result = BDN_XOR_DEC_DEC(manager, f, mux_dec, T);
            __DSD_RecursiveDeref(manager, mux_dec);        
            return result;


        
        }
    }

    if((GET_TYPE(TReg)) == DSD_XOR)
    {	
        if(node_exists(TReg->actual_list, E))
        {
#ifdef DEBUG_PRINT
            printf("-------Section 2.5.1----------\n");
#endif

            temp = (ActualNode*) FixHeapMalloc(actual_malloc_ptr);
            temp->next = NULL;
            temp->decomposition = E;


            Tresidue = list_residue(manager->Ddmanager_analogue, TReg->actual_list, temp, size);
            Tresidue_dec = BDN_BDD_XOR_RESIDUE(manager, Tresidue);


            free(size);
            FixHeapFree(actual_malloc_ptr, temp);

            if(!DSD_IsComplement(T))
            {
                mux_dec = BDN_BDD_MUX_VAR_DEC_DEC(manager, top_func, DSD_Not(manager->one), Tresidue_dec);
                __DSD_RecursiveDeref(manager, Tresidue_dec);

                result =  BDN_XOR_DEC_DEC(manager, f, mux_dec, E);
                __DSD_RecursiveDeref(manager, mux_dec);        
                return result;


            
            }
            else
            {
                mux_dec = BDN_BDD_MUX_VAR_DEC_DEC(manager, top_func, manager->one, Tresidue_dec);
                __DSD_RecursiveDeref(manager, Tresidue_dec);

                result = BDN_NXOR_DEC_DEC(manager, f, mux_dec, E);
                __DSD_RecursiveDeref(manager, mux_dec);        
                return result;            
            }	


        }

#ifdef DEBUG_PRINT
        printf("-------Section 2.5.2----------\n");
#endif
    }

    free(size);
    return NULL;

}





