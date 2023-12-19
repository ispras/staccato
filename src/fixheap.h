/*
 * This file was written by Todd M. Austin, University of Wisconsin - Madison
 * Department of Computer Sciences
 *
 * Copyright (C) 1992, 1993, 1994, 1995 by Todd M. Austin
 *
 * This source file is distributed in the hope that it will be useful,
 * but without any warranty.  No author or distributor accepts
 * responsibility to anyone for the consequences of using it or for
 * whether it serves any particular purpose or works at all.
 *
 * Everyone is granted permission to copy, modify and redistribute
 * this source file under the following conditions:
 *
 *    Permission is granted to anyone to make or distribute copies
 *    of this source code, either as received or modified, in any
 *    medium, provided that all copyright notices, permission and
 *    nonwarranty notices are preserved, and that the distributor
 *    grants the recipient permission for further redistribution as
 *    permitted by this document.
 *
 *    Permission is granted to distribute this file in compiled
 *    or executable form under the same conditions applying for
 *    source code, provided that either
 *    A. it is accompanied by the corresponding machine-readable
 *       source code, or
 *    B. it is accompanied by a written offer, with no time limit,
 *       to give anyone a machine-readable copy of the corresponding
 *       source code in return for reimbursement of the cost of
 *       distribution.  This written offer must permit verbatim
 *       duplication by anyone.
 *    C. it is distributed by someone who received only the
 *       executable form, and is accompanied by a copy of the
 *       written offer of source code which s/he received along with it.
 *
 * In other words, you are welcome to use, share and improve this
 * source file.  You are forbidden to forbid anyone else to use, share
 * and improve what you give them.
 *
 * INTERNET: austin@cs.wisc.edu
 * US Mail: 1210 W. Dayton Street, Madison, WI 53706
 *
 * $Id: fixheap.h,v 1.1.1.1 2004/07/03 17:58:17 qedq Exp $
*/

#ifndef FIXHEAP_H
#define FIXHEAP_H

#include <stdio.h>
#define MALLOC_OVERHEAD 64

/*!
  fixed size node heap manager:
 
  most fast Unix malloc packages allocate blocks in only sizes that are
  a power of 2 bytes; for a program that uses a LOT of small structures
  (like TETRA) this greatly increase memory usage; for example, a
  38 byte structure will typically be allocated in a 64 byte block; the
  fixed node size heap manager allocates blocks in a large array, and then
  maintains a free list for this allocation, thus the system malloc overhead
  is only bore for a block of allocations, rather than for each; the other
  advantage of a fixed size allocator is better memory reference locality -
  since the overhead is moved away from the allocated nodes
 
  do not use this for infrequently used nodes; to compute the nodesPerAlloc,
  use a variant of the following equation
 
  	nodesPerAlloc = floor((64k - sys_overhead - FixHeap_overhead)/nodeSize)
 
  where:
 	sys_overhead: is the system malloc overhead (always less 64 bytes on
 	  the DECstation's default malloc package), make this conservatively
 	  large!
 	FixHeap_overhead: is the fixed head size manager's overhead, which
 	  is always less than 32 bytes
 
  Increase or decrease the mass allocation size to match the usage of the
  particular node type.
*/

typedef struct _FixHeap {
    char *name;
    unsigned long nodeSize; /* note: nodeSize>=sizeof(void *), for freeList */
    int verbose;
#ifdef FIXED_SIZE_HEAP
    unsigned long nodesPerAlloc;
    void *freeList;
    void *heapAlloc;    /* list of allocated blocks */
#endif /* FIXED_SIZE_HEAP */
} FixHeap, *FixHeapPtr;

FixHeapPtr FixHeapCreate(char *name, unsigned long nodeSize, 
        unsigned long nodesPerAlloc, int verbose);
void FixHeapRelease(FixHeapPtr heap);
void *FixHeapMalloc(FixHeapPtr heap);
void FixHeapFree(FixHeapPtr heap, void *node);

/*!
  use this to compute reasonable malloc chunk sizes 
*/
#define FIXHEAP_OVERHEAD	32

/*!
  compute the number of nodes to allocate per heap expand, this is
  specified as the size (in kBytes) of each chunk allocated, for example,
  COMPUTE_ALLOC_COUNT(32, FooRec), will make the fix heap manager allocate
  as many FooRec's that fit into 32k each time it expands the heap, the
  fix heap manager knows about the overheads incurred by malloc operations
  (from machine.h) and fix heap management (FIXHEAP_OVERHEAD)
*/
#define COMPUTE_ALLOC_COUNT(kBytes,Type) \
(((((kBytes)*1024)-MALLOC_OVERHEAD)-FIXHEAP_OVERHEAD)/sizeof(Type))

#endif /* FIXHEAP_H */
