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


#include <stdio.h>
#include "DSDManager.h"

FixHeapPtr dsd_malloc_ptr;
FixHeapPtr actual_malloc_ptr;

/*table maintenance stuff*/
DSDNode **Create_DSD_Table(int size)
{
    int i;
    DSDNode **iter;

    assert((iter = (DSDNode **) malloc((size * sizeof(DSDNode*)))) != NULL);


    for(i = 0; i < size; i++)
    {
        iter[i] = NULL;
    }

    return iter;
}


void delete_actual_list(DSDManager *manager, DSDNode *dsd)
{
    int i;
    ActualNode *iter;
    ActualNode *iter2;

    iter = DSD_Regular(dsd)->actual_list;

    while(iter)
    {
        iter2 = iter;
        iter = iter2->next;

        FixHeapFree(actual_malloc_ptr, iter2);

    }
}

void delete_support(DSDNode *dsd)
{
    free(DSD_Regular(dsd)->support);
}



void Destroy_DSD_Table(DSDManager *manager, DSDNode **DSDTable, int size)
{
    int i;
    DSDNode * iter;
    DSDNode * iter_copy;

    manager->dead_nodes_threshold = MAX_NOT_ALLOWED;

    for(i = 0; i < size; i++)
    {
        iter = DSDTable[i];

        while(iter)
        {
            assert(!DSD_IsComplement(iter));

            iter_copy = iter->next;

            if(GET_TYPE(iter) != VAR) Cudd_RecursiveDeref(manager->Ddmanager_analogue, iter->bdd_analogue);
#ifndef DISABLE_SBDD            
            if(GET_TYPE(iter) != VAR) Cudd_RecursiveDeref(manager->Ddmanager_analogue, iter->symbolic_kernel);
#endif
            delete_actual_list(manager, iter);
            //delete_support(iter);

            FixHeapFree(dsd_malloc_ptr, iter);
            iter = iter_copy;
        }
    }


    FixHeapRelease(dsd_malloc_ptr);
    FixHeapRelease(actual_malloc_ptr);
    free(DSDTable);

}

void purge_triggered_stat_update(DSDManager *manager)
{
    return;
    /*
    
    if(manager->total_actualsize > manager->max_actualsize)
    {
        manager->max_actualsize = manager->total_actualsize;
    }

    manager->current_average_actualsize = ((double) manager->total_actualsize) / (manager->num_DSD_nodes);    

    if(manager->current_average_actualsize > manager->max_average_actualsize)
    {
        manager->max_average_actualsize = manager->current_average_actualsize;
    }


    manager->current_memory_used = manager->total_actualsize * sizeof(ActualNode) + manager->num_DSD_nodes * sizeof(DSDNode) + manager->support_size;

    if(manager->current_memory_used > manager->max_memory_used)
    {
        manager->max_memory_used = manager->current_memory_used;
    }

    if(manager->num_DSD_nodes > manager->max_DSD_nodes)
    {
        manager->max_DSD_nodes = manager->num_DSD_nodes;
    }

    if(manager->support_size > manager->max_support_size)
    {
        manager->max_support_size = manager->support_size;
    }

    */

}    


/*possible short-circuit when the correct number
  of elements are freed*/
void Purge_Derefs(DSDManager *manager)
{
    int i;
    DSDNode * iter;
    DSDNode * iter_copy;
    DSDNode * previous;
    int support_size;
    int num;

    num = 0;
    
    manager->garbage_cleans++;
   
    manager->dead_nodes_threshold += manager->dead_nodes_threshold;
    
    assert(manager->dead_nodes_threshold != MAX_NOT_ALLOWED);

    purge_triggered_stat_update(manager);

    for(i = 0; i < manager->DSD_unique_table_size; i++)
    {
        iter = manager->DSD_unique_table[i];
        previous = iter;

        while(iter != NULL)
        {
            assert(!DSD_IsComplement(iter));

            if(GET_REF(iter) == 0 && (GET_TYPE(iter) != VAR))
            {
                num++;
                assert(iter->support);
                
                iter_copy = iter->next;

                if(iter == manager->DSD_unique_table[i])
                {
                    manager->DSD_unique_table[i] = iter_copy;
                }
                else
                {
                    previous->next = iter_copy;
                }

                /*probably will be unnecessary in the future*/
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, iter->bdd_analogue);
#ifndef DISABLE_SBDD            
                Cudd_RecursiveDeref(manager->Ddmanager_analogue, iter->symbolic_kernel);
#endif

                /*
                if((DSD_Regular(iter))->support)
                {
                    support_size = (GET_CAN((DSD_Regular(iter))) / (sizeof(int) * 8)) + 1;
                    manager->support_size -= (support_size*sizeof(int));
                }
                */ 
                
//                manager->total_actualsize -= INPUT_SIZE((DSD_Regular(iter)));
                

                delete_actual_list(manager, iter);
                /*delete_support(iter);*/

                FixHeapFree(dsd_malloc_ptr, iter);

                iter = iter_copy;

            }
            else
            {
                previous = iter;
                iter = iter->next;
            }
        }
    }



    manager->num_DSD_nodes = (manager->num_DSD_nodes) -	num;
    manager->current_average_actualsize = ((double)(manager->total_actualsize))/(manager->num_DSD_nodes);
    manager->current_memory_used = manager->total_actualsize * sizeof(ActualNode) + manager->num_DSD_nodes * sizeof(DSDNode) + manager->support_size;

    manager->dead_nodes_current = 0;

}

void __DSD_RecursiveDeref(DSDManager *manager, DSDNode * dsd)
{
    recursive_deref(manager, dsd);
    if(manager->dead_nodes_current > manager->dead_nodes_threshold)
    {
        Purge_Derefs(manager);   	
    }
}

void recursive_deref(DSDManager *manager, DSDNode *dsd2)
{
    int i;
    ActualNode * temp_actual;
    DSDNode *dsd;
    
    if(!dsd2) return;

    dsd = DSD_Regular(dsd2);

    temp_actual = dsd->actual_list;

    if(GET_REF(dsd) < 1)
    {
        return;
    }


    assert(GET_REF(dsd) > 0);

    if(GET_REF(dsd) != SATURATION && GET_TYPE(DSD_Regular(dsd)) != VAR)
    {
        dsd->topvar_refsize = ((dsd->topvar_refsize) & 0xffff0000) | ((dsd->topvar_refsize & 0x0000ffff) - 1);
    }

    if(GET_REF(dsd) == 0 && GET_TYPE(DSD_Regular(dsd)) != VAR)
    {
        dsd->support = (int*) 1;
        
        manager->dead_nodes_current = manager->dead_nodes_current + sizeof(DSDNode) + sizeof(ActualNode)*INPUT_SIZE(dsd);
        while(temp_actual)
        {
            recursive_deref(manager, temp_actual->decomposition);
            temp_actual = temp_actual->next;
        }
    }

}

/*refencing will only be called by the decomposition algorithm
  make sure that the ref is not double incremented*/
void __DSD_Ref(DSDManager *manager, DSDNode * dsd)
{
    ActualNode *iter;
    
    if(GET_REF(DSD_Regular(dsd)) != SATURATION)	
        (DSD_Regular(dsd)->topvar_refsize)++;
    

    if(GET_REF(DSD_Regular(dsd)) == 1 && DSD_Regular(dsd)->support)
    {
        DSD_Regular(dsd)->support = NULL;
        
        manager->dead_nodes_current = manager->dead_nodes_current - sizeof(DSDNode) - sizeof(ActualNode)*INPUT_SIZE(DSD_Regular(dsd));
        
        iter = (ActualNode *) DSD_Regular(DSD_Regular(dsd)->actual_list);
        
        while(iter)
        {
            __DSD_Ref(manager, iter->decomposition);
            iter = (ActualNode *) DSD_Regular(iter->next);
        }
    }
        
}


DSDNode * find_DSD_node(DSDManager *manager, DdNode *bdd)
{
    unsigned int index;
    int table_size;

    DSDNode *result;

    if(Cudd_Regular(bdd) == Cudd_ReadOne(manager->Ddmanager_analogue))
    {
        if(Cudd_IsComplement(bdd)) return DSD_Not(manager->one);
        return manager->one;
    }

    index = (long int) Cudd_Regular(bdd);
    table_size = manager->DSD_unique_table_size;

    result = manager->DSD_unique_table[((index>>4)%table_size)];

    
    while(result)
    {
        assert(!DSD_IsComplement(result));
        if(Cudd_Not(result->bdd_analogue) == bdd)
        {
            if(DSD_Regular(result)->support)
            {
                __DSD_Ref(manager, result);
                DSD_Regular(result)->topvar_refsize = ((DSD_Regular(result)->topvar_refsize) & 0xffff0000) | ((DSD_Regular(result)->topvar_refsize & 0x0000ffff) - 1);
            }
            
            return DSD_Not(result);
        }
        else if(result->bdd_analogue == bdd)
        {
            if(DSD_Regular(result)->support)
            {
                __DSD_Ref(manager, result);
                DSD_Regular(result)->topvar_refsize = ((DSD_Regular(result)->topvar_refsize) & 0xffff0000) | ((DSD_Regular(result)->topvar_refsize & 0x0000ffff) - 1);
            }
            return result;
        }

        result = result->next;
    }

    
    
    return result;

}

DSDNode *create_DSD_node(DSDManager *manager, DdNode *bdd)
{
    unsigned int index;
    int table_size;
    int hash_index;

    DSDNode *dsd_node;
    DSDNode *temp_node;


    assert(Cudd_Regular(bdd) != manager->one->bdd_analogue); 

    index = (long int) Cudd_Regular(bdd);
    table_size = manager->DSD_unique_table_size;

    hash_index = (index>>4)%table_size;

    manager->num_DSD_nodes++;

    dsd_node = (DSDNode*) FixHeapMalloc(dsd_malloc_ptr);

    dsd_node->next = NULL;
    dsd_node->bdd_analogue = bdd;
#ifndef DISABLE_SBDD
    dsd_node->symbolic_kernel = NULL;
#else
    dsd_node->symbolic_kernel = bdd;
#endif
    dsd_node->actual_list = NULL;
    dsd_node->topvar_refsize = 0;
    dsd_node->support = NULL;
    dsd_node->parent = NULL;


    Cudd_Ref(dsd_node->bdd_analogue);

    dsd_node->type_actualsize = 0;

    /*calling function is responsible for setting the type*/

    temp_node = manager->DSD_unique_table[hash_index];


    dsd_node->next = temp_node;

    manager->DSD_unique_table[hash_index] = dsd_node;
    assert(!DSD_IsComplement(temp_node));
    assert(!DSD_IsComplement(dsd_node));

    return dsd_node;

}

void __DSD_Quit(DSDManager *manager)
{
    /*destroy unique table*/
    Destroy_DSD_Table(manager, manager->DSD_unique_table, manager->DSD_unique_table_size);

    /*destroy other tables--recommendation cache, supports, etc*/
    //FixHeapFree(dsd_malloc_ptr, manager->one);

    /*delete manager*/
    free(manager);

   /* 
    map<unsigned int, int>::iterator iter;

    for(iter = REF_COUNT.begin(); iter != REF_COUNT.end(); iter++)
    {
      printf("BDD: %x   Reference Count: %d    Final Count: %d\n", (*iter).first, (*iter).second, ((DdNode*)((*iter).first))->ref);
    }
    */
}



DSDManager* __DSD_Init(DdManager* manager, int_32 recommendation_size)
{
    DSDManager *dmanager;
    int unique_size;  /*allow to be parameterized in the future*/
    int nodesPerAlloc;

    DSDNode* one;
    /*eventually support dynamic initialization--currently
      only one manager is allowed*/


    nodesPerAlloc = COMPUTE_ALLOC_COUNT(16,DSDNode);
    dsd_malloc_ptr = FixHeapCreate((char *)"dsd_malloc", sizeof(DSDNode), nodesPerAlloc, 0);

    nodesPerAlloc = COMPUTE_ALLOC_COUNT(16,DSDNode);
    actual_malloc_ptr = FixHeapCreate((char *)"actual_malloc", sizeof(ActualNode), nodesPerAlloc, 0);


    one = (DSDNode*) FixHeapMalloc(dsd_malloc_ptr);

    one->bdd_analogue = Cudd_ReadOne(manager);
    one->symbolic_kernel = NULL;
    one->type_actualsize = 0;
    one->actual_list = NULL;
    one->support = NULL;
    one->parent = NULL;
    one->next = NULL;
    one->topvar_refsize = 0;

    dmanager = (struct DSDManager *) malloc(sizeof(struct DSDManager));

    dmanager->one = one;
    dmanager->notdisjoint = 0;
    dmanager->num_disjoint = 0;

    
    dmanager->decomposed_outputs = 0;
    dmanager->num_outputs = 0;
    dmanager->num_blocks = 0;

    dmanager->num_entered = 0;
    dmanager->num_primes = 0;
    dmanager->num_commons = 0;
    dmanager->num_newdecomp = 0;
    dmanager->support_size = 0;
    dmanager->num_unique_blocks = 0;
    dmanager->num_unique_symbolic_blocks = 0;
    dmanager->theoretical_DSD_consumption = 0;
    dmanager->theoretical_Actual_consumption = 0;
    dmanager->theoretical_memory_consumption = 0;


    /*
       create unique table and initial variables, mem
       create symbolic manager
       create and initialize array formal  
       create and initialize added bdds




       create other caches

     */

    unique_size = 24593; /*support ability for this to double*/
    
    dmanager->DSD_unique_table = Create_DSD_Table(unique_size);
    dmanager->DSD_unique_table_size = unique_size;
    dmanager->num_DSD_nodes = 0;
    
    if(recommendation_size < 1)
    {
      dmanager->dead_nodes_threshold = RECOMMENDATION_DEFAULT;
    }
    else
    {  
      dmanager->dead_nodes_threshold = recommendation_size;  /*this should dynamically change*/
    }

    dmanager->dead_nodes_current = 0;

    dmanager->garbage_cleans = 0;
    
    /*dmanager->recommendation_cache = initialize_recommendation_cache(recommendation_size);*/
    
    
    dmanager->max_DSD_nodes = 0;

    dmanager->max_support_size = 0;

    dmanager->max_average_actualsize = 0;
    dmanager->max_actualsize = 0;
    dmanager->total_actualsize = 0;
    dmanager->current_average_actualsize = 0;

    dmanager->current_memory_used = 0;
    dmanager->max_memory_used = 0; 

    dmanager->Ddmanager_analogue = manager;
    
    dmanager->snode_size = 0;
    dmanager->onode_size = 0;
    dmanager->snode_counter = 0;
    dmanager->onode_counter = 0;
    
    dmanager->node_size = 0;
    dmanager->node_counter = 0;
    return dmanager;
}

