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


#ifndef _DSD
#define _DSD

#ifdef DISABLE_SBDD
  #define DISABLE_GC
  #define DISABLE_SM
#endif

#include "util.h"
#include "cudd.h"
#include "cuddInt.h"
#include "fixheap.h"
#include <assert.h>
#include <stdio.h>


/*!
  Different DSD types possible
*/
#define VAR 0
#define PRIME 1
#define OR  2
#define XOR 3

/*!
  Internal defines
*/
#define MAX_NOT_ALLOWED ~0
#define TYPE_MASK 127<<25
#define SIZE_MASK 255<<24
#define MARK_MASK 1<<24
#define SATURATION 16384
#define CAN_MASK ((1<<16) - 1)<<16 



typedef struct DSDNode DSDNode;
typedef unsigned int int_32;
typedef unsigned short int_16;
typedef struct ActualNode ActualNode;
typedef struct DSDManager DSDManager;
typedef struct SupportList SupportList;



extern FixHeapPtr dsd_malloc_ptr;
extern FixHeapPtr actual_malloc_ptr;

/*! \struct
  Linked list node that contains a pointer to a decomposition
  and a pointer to the next AcualNode in the actuals list.
*/
struct ActualNode{
  /*!
    Points to a decomposition
  */
  DSDNode * decomposition;
  /*!
    Points to the next ActualNode in the list--NULL
    if it is the last entry
  */
  ActualNode * next;
};

/*! \struct
  This structure is no longer used in STACCATO.  Previously,
  it was used to explicitly store the support of each DSDNode
*/
struct SupportList{
  int_32 support_var;
  SupportList *next;
};



/*! \struct
  The fundamental structure of STACCATO.  This structure
  is comparablef to the DdNode in CUDD.  It stores information
  on decomposition type and the actuals list for the decomposition.
  It also stores internal reference infomation along with pointers
  to the coresponding BDD and symbolic kernel
*/
struct DSDNode{
  /*!
    This variable encodes information for decomposition type
    and the size of actuals list
  */
  int_32 type_actualsize;
  
  /*!
    This variable contains information on the decomposition
    lowest BDD variable along with the reference count
    for this DSD
  */ 
  int_32 topvar_refsize;

  /*!
    This a linked list of ActualNodes, i.e., this is the actuals list
  */
  ActualNode * actual_list;
  
  /*!
    Pointer to the symbolic kernel function for this DSD
  */
  DdNode * symbolic_kernel;
  
  /*!
    Pointer to the BDD needed to build the decomposition and used
    for indexing purposes
  */
  DdNode * bdd_analogue;
  
  /*!
    Pointer to the next DSDNode in the unique table that holds all the DSDNodes
  */ 
  DSDNode * next;
  
  DSDNode *parent;
  int *support;
};


/*! \struct
  This structure is an analogue to the DdManager.  It contains
  the location of the corresponding DdManager along with cache
  and initialization information.
*/
struct DSDManager{
  /*!
    Pointer to the corresponding DdManager
  */
  DdManager * Ddmanager_analogue;

  /*!
    Represents the constant 1 node
  */
  DSDNode *one;	

  /*!
    Number of times DSD_Create is called
  */
  int_32 num_outputs;
  
  /*!
    Number of DSD produced by DSD_Create that
    have some decomposability
  */
  int_32 decomposed_outputs;
  
  /*!
    Number of blocks across all the DSDs currently
    referenced in the DSDManager after calling
    the function update_blocks.
  */
  int_32 num_blocks;
  
  int_32 num_unique_blocks;
  int_32 num_unique_symbolic_blocks;
  int_32 theoretical_DSD_consumption; 
  int_32 theoretical_Actual_consumption;
  int_32 theoretical_memory_consumption;
  int_32 garbage_cleans;

  /*!
    Unique table for DSDs
  */
  DSDNode ** DSD_unique_table;
  
  /*!
    Current size of the DSD unique table
  */
  int_32 DSD_unique_table_size;

  int_32 notdisjoint;

  int_32 num_DSD_nodes;
  int_32 max_DSD_nodes;

  int_32 support_size;
  int_32 max_support_size;

  double max_average_actualsize;
  double current_average_actualsize;
  int_32 max_actualsize;
  int_32 total_actualsize;

  int_32 max_memory_used; 
  int_32 current_memory_used;  


  /*!
    The number of dereferenced DSDNodes
    required to trigger garbage collection.  This
    number increase with every garbage collection.
  */
  int_32 dead_nodes_threshold; 

  /*!
    The current number of unreferenced DSDNodes.
  */
  int_32 dead_nodes_current; 


  int_32 num_disjoint;
  int_32 num_commons;
  int_32 num_entered;
  int_32 num_primes;
  int_32 num_newdecomp;

  int_32 num_snodes;
  int_32 num_onodes;
  DdNode **snodes_array;
  DdNode **onodes_array;
  int_32 snode_size;
  int_32 onode_size;
  int_32 snode_counter;
  int_32 onode_counter;

  int_32 num_nodes;
  DdNode **nodes_array;
  int_32 node_size;
  int_32 node_counter;

};


#define DSD_Regular(x) ((DSDNode *)((long int)(x) & ~01))
#define DSD_Not(x) ((DSDNode *) ((long int)(x) ^ 01))
#define DSD_Complement(x) ((DSDNode *) ((long int)(x) | 01))
#define DSD_IsComplement(x) ((int) ((long int) (x) & 01))


#define GET_TYPE(x) (((x)->type_actualsize & TYPE_MASK)>>25)
#define SET_TYPE(x, type) ((x)->type_actualsize = (((x)->type_actualsize & ~(TYPE_MASK)) | (type << 25)))
#define INPUT_SIZE(x) ((x)->type_actualsize & ~(SIZE_MASK))
#define SET_SIZE(x, size) ((x)->type_actualsize = (((x)->type_actualsize & (SIZE_MASK)) | (size)))
#define mark(x) ((x)->type_actualsize = (((x)->type_actualsize | (MARK_MASK))))
#define unmark(x) ((x)->type_actualsize = (((x)->type_actualsize & ~(MARK_MASK))))

#define marked(x) (((x)->type_actualsize & (MARK_MASK))>>24)
#define SET_CAN(x, value) ((x)->topvar_refsize = (((x)->topvar_refsize & ~(CAN_MASK)) | (value << 16)))
#define GET_CAN(x) (((x)->topvar_refsize & CAN_MASK)>>16)
#define GET_REF(x) ((x)->topvar_refsize & ~(CAN_MASK))

#endif
