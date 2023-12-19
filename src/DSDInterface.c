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


#include "DSDInterface.h"





void print_stats(DSDManager* manager)
{
  /*	printtf("Total Number of Outputs: %d\n", manager->num_outputs);
        printf("Total Decomposed Outputs: %d\n", manager->decomposed_outputs);
        printf("Total Number of Blocks: %d\n", manager->num_blocks);
        printf("Total Number of Non-trivial Decompositions: %d\n", manager->num_entered);
        printf("Times Prime Algorithm Used: %d\n", manager->num_primes);
        printf("Times Common Algorithm Used: %d\n", manager->num_commons);
        printf("Times New Decomposition Algorithm used: %d\n", manager->num_newdecomp);

        printf("Number of DSD Nodes created: %d\n\n\n", manager->num_DSD_nodes);
   */


  purge_triggered_stat_update(manager);

  /*
     printf("%d & ", manager->num_outputs);
     printf("%d & ", manager->decomposed_outputs);
     printf("%d & ", manager->num_blocks);
     printf("%d & ", manager->num_entered - (manager->num_primes + manager->num_newdecomp));
     printf("%.1f & ", (((float) manager->num_primes)/manager->num_entered) * 100);
     printf("%.1f \\\\\n", (((float) manager->num_newdecomp)/manager->num_entered) * 100);
   */

  printf("Number of outputs: %d\n", manager->num_outputs);
  printf("Number of Decomposed outputs: %d\n", manager->decomposed_outputs);
  printf("Number of Blocks Used: %d\n", manager->num_blocks);

  //printf("Number of Unique Blocks Used: %d\n", manager->num_unique_blocks);

  //printf("Number of Unique Symbolic Blocks Used: %d\n", manager->num_unique_symbolic_blocks);

  //printf("Percentage of symbolic reduction: %.1f \n", ((float) (manager->num_unique_blocks - manager->num_unique_symbolic_blocks))/((float) manager->num_unique_blocks)*100); 

  //printf("Number of times algorithm entered: %d \n", manager->num_entered);
  printf("Percentage of Prime Decompositions: %.1f\n", (((float) manager->num_primes)/manager->num_entered) * 100);
  printf("Percentage of New Decompositions %.1f \n", (((float) manager->num_newdecomp)/manager->num_entered) * 100);
  printf("Percentage of 4.1 Decompositions %.1f \n", (((float) manager->num_disjoint)/manager->num_entered) * 100);
  //printf("Times New Decomposition Algorithm used: %d\n\n\n", manager->num_newdecomp);

  //printf("Speculative Final DSD Nodes created: %d\n", manager->num_DSD_nodes);

  //printf("Final DSD Nodes created: %d\n", manager->max_DSD_nodes);
  //printf("Final Actual Nodes created: %d\n", manager->max_actualsize);
  //printf("Final Average Actual List size: %.1f\n\n", ((float) manager->max_actualsize)/(manager->max_DSD_nodes));


  //printf("Number of times garbage collected: %d \n", manager->garbage_cleans - 1);
  //printf("Size of the free node cache: %d \n", manager->dead_nodes_threshold);


  //printf("Final DSD Memory Consumed: %d BYTES\n", (manager->max_DSD_nodes)*sizeof(DSDNode));
  //printf("Final Actual Node Memory Consumed: %d BYTES\n", (manager->max_actualsize)*sizeof(ActualNode));


  manager->max_memory_used = (manager->max_DSD_nodes)*sizeof(DSDNode) + (manager->max_actualsize)*sizeof(ActualNode);

  /*printf("Max Support Memory Consumed: %d BYTES\n", manager->max_support_size);*/
  /*printf("Max Total Memory Consumed: %d BYTES\n\n\n", (manager->max_memory_used));*/

  /*
     printf("Current DSD Nodes created: %d\n", manager->num_DSD_nodes);
     printf("Current Actual Nodes created: %d\n", manager->total_actualsize);
     printf("Current Average Actual List size: %f\n\n", manager->current_average_actualsize);


     printf("Current DSD Memory Consumed: %d BYTES\n", (manager->num_DSD_nodes)*sizeof(DSDNode));
     printf("Current Actual Node Memory Consumed: %d BYTES\n", (manager->total_actualsize)*sizeof(ActualNode));
     printf("Current Support Memory Consumed: %d BYTES\n", manager->support_size);
     printf("Current Total Memory Consumed: %d BYTES\n", (manager->current_memory_used));
   */

}


/*DSD manager functions--call in sync with DD*/

DSDManager* DSD_Init(DdManager* manager, int_32 recommendation_size)
{
  return __DSD_Init(manager, recommendation_size);
}


void update_blocks(DSDManager *manager, DSDNode *result)
{
  ActualNode *iter;

  if(DSD_Regular(result) == manager->one)
  {
    return;
  }

  manager->onode_size++;
  manager->snode_size++;

  if(GET_TYPE(DSD_Regular(result)) == VAR)
  {
    return;
  }

  /*    if((GET_TYPE(DSD_Regular(result)) == OR) || (GET_TYPE(DSD_Regular(result)) == XOR))
        {
        manager->onode_size--;
        manager->snode_size--;
        }    
   */    
  if(GET_TYPE(DSD_Regular(result)) == OR || GET_TYPE(DSD_Regular(result)) == XOR)
  {
    manager->num_blocks += (INPUT_SIZE(DSD_Regular(result)) - 1);
  } 
  else
  {
    manager->num_blocks++;
  }

  iter = DSD_Regular(result)->actual_list;

  while(iter)
  {
    update_blocks(manager, iter->decomposition);
    iter = iter->next;
  }

}

void load_arrays(DSDManager *manager, DSDNode *result)
{
  ActualNode *iter;

  if(DSD_Regular(result) == manager->one)
  {
    return;
  }



  /*    if(!((GET_TYPE(DSD_Regular(result)) == OR) || (GET_TYPE(DSD_Regular(result)) == XOR)))
        {*/
  manager->snodes_array[manager->snode_counter++] = DSD_Regular(result)->symbolic_kernel;

  manager->onodes_array[manager->onode_counter++] = DSD_Regular(result)->bdd_analogue;

  manager->nodes_array[manager->node_counter++] = DSD_Regular(result)->bdd_analogue;
  manager->nodes_array[manager->node_counter++] = DSD_Regular(result)->symbolic_kernel;
  /*    }*/

  iter = DSD_Regular(result)->actual_list;

  while(iter)
  {
    load_arrays(manager, iter->decomposition);
    iter = iter->next;
  }
}

void validate_deref(DSDManager *manager)
{
  DSDNode *iter;
  int i;

  for(i = 0; i < manager->DSD_unique_table_size; i++)
  {
    iter = manager->DSD_unique_table[i];

    while(iter != NULL)
    {
      if(!marked(iter)) 
        assert(0); 
      iter = iter->next;
    }
  }

}

void count_unique(DSDManager *manager, DSDNode *result)
{
  ActualNode *iter;
  int symbolic;
  DSDNode *temp;

  symbolic = 0;

  if(DSD_Regular(result) == manager->one)
  {
    return;
  }

  if(marked(DSD_Regular(result)))
  {
    return;
  }

  mark(DSD_Regular(result));

  //manager->max_DSD_nodes++; 

  if(GET_TYPE(DSD_Regular(result)) == VAR)
  {
    return;
  }

  manager->total_actualsize+= INPUT_SIZE(DSD_Regular(result));
#if 0
  temp = find_DSD_node(manager, DSD_Regular(result)->symbolic_kernel);

  if(!temp)
  {
    symbolic = 1;
    temp = create_DSD_node(manager, DSD_Regular(result)->symbolic_kernel);    
    manager->num_DSD_nodes--;
    Cudd_Ref(DSD_Regular(result)->symbolic_kernel);
    DSD_Regular(temp)->support = (int*) 1;
    DSD_Regular(temp)->symbolic_kernel = DSD_Regular(result)->symbolic_kernel;
  } 
  else if(!(DSD_Regular(temp)->support))
  {
    /* if(marked(DSD_Regular(temp)))
       {
       symbolic = 1;
       }*/ 
    symbolic = 1; 
    DSD_Regular(temp)->support = (int*) 1;
  }



  manager->max_actualsize += INPUT_SIZE(DSD_Regular(result));

  if(GET_TYPE(DSD_Regular(result)) == OR || GET_TYPE(DSD_Regular(result)) == XOR)
  {
    manager->num_unique_blocks += (INPUT_SIZE(DSD_Regular(result)) - 1);

    if(symbolic)
    {
      DSD_Regular(result)->support = (int*) 1;
      manager->num_unique_symbolic_blocks += (INPUT_SIZE(DSD_Regular(result)) - 1);
    } 

  } 
  else
  {
    manager->num_unique_blocks++;

    if(symbolic)
    {
      DSD_Regular(result)->support = (int*) 1;
      manager->num_unique_symbolic_blocks++;
    } 
  }

#endif
  iter = DSD_Regular(result)->actual_list;

  while(iter)
  {
    count_unique(manager, iter->decomposition);
    iter = iter->next;
  }

}



void DSD_Quit(DSDManager *manager)
{
  __DSD_Quit(manager);
}

/*analogues to the ref and deref functions in cudd*/
DSDNode* DSD_Create(DSDManager* manager, DdNode* f)
{
  DSDNode *result;


  manager->num_outputs++;

  result = __DSD_Create(manager, f);

  if(GET_TYPE(DSD_Regular(result)) == PRIME) /*|| GET_TYPE(DSD_Regular(result)) == OR || GET_TYPE(DSD_Regular(result)) == XOR)*/
  {
#ifndef DISABLE_SBDD
    if((DSD_Regular(result)->bdd_analogue != DSD_Regular(result)->symbolic_kernel))
    {
      manager->decomposed_outputs++;
    }
    else
    {
      //printf("No decomposition\n");
    }
#endif
  }
  else if(DSD_Regular(result) != manager->one)
  {
    manager->decomposed_outputs++;
  }


  //update_blocks(manager, result);


  return result;

}




void DSD_RecursiveDeref(DSDManager *manager, DSDNode * dsd_node)
{
  __DSD_RecursiveDeref(manager, dsd_node);
}

void DSD_Deref(DSDManager *manager, DSDNode * dsd_node)
{
  __DSD_RecursiveDeref(manager, dsd_node);
}

void DSD_Ref(DSDManager *manager, DSDNode * dsd_node)
{
  __DSD_Ref(manager, dsd_node);
}




int_32 Get_Type(DSDNode* dsd_node)
{
  return GET_TYPE(dsd_node);
}






int_32 get_num_actuals(DSDNode* dsd_node)
{
  return __Get_Input_Count(dsd_node);
}

/*
DSDNode* Get_First_Input(DSDNode* dsd_node)
{
  return __Get_First_Input(dsd_node);
}
*/

DSDNode *Get_X_Input(DSDNode* dsd_node, int index)
{
  int i;
  DSDNode *dsd_node_reg;
  ActualNode *actual_iter;
  if(index < 0) return NULL;
  
  dsd_node_reg = DSD_Regular(dsd_node);
  
  if(GET_TYPE(dsd_node_reg) == VAR) return NULL;
  if(index >= INPUT_SIZE(dsd_node_reg)) return NULL;

  i = 0;
  actual_iter = dsd_node_reg->actual_list; 
  
  
  while(i < index)
  {
    actual_iter = actual_iter->next;
    i++;
  }
 
  return actual_iter->decomposition;
}
 
void mark_decomposition(DSDNode *dsd_node)
{
  mark(DSD_Regular(dsd_node));
}


void unmark_decomposition(DSDNode *dsd_node)
{
  unmark(DSD_Regular(dsd_node));
}

int is_marked(DSDNode *dsd_node)
{
  return marked(DSD_Regular(dsd_node));
}

void Decomposition_Print(DSDNode * dsd_node)
{
  __Decomposition_Print(dsd_node);
}

void Recursive_Decomposition_Print(DSDNode * dsd_node)
{
  __Recursive_Decomposition_Print(dsd_node);
}

DdNode * get_symbolic_kernel(DSDNode* dsd_node)
{
  return __Get_Symbolic_Decomposition(dsd_node);
}

DdNode * get_bdd(DSDNode* dsd_node)
{
  return __Get_BDD(dsd_node);
}


