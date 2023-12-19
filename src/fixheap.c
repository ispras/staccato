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
 * $Id: fixheap.c,v 1.1.1.1 2004/07/03 17:58:17 qedq Exp $
 */


#include <stdlib.h>

/*
 *        various compilation options
 */

/* define INLINE to 'static' if your system cannot handle 'inline', all
   C++ compilers and the GNU gcc compiler support the keyword "inline",
   this has been moved to the makefile
 */
#ifndef INLINE
#ifdef __GNUC__
#define INLINE static inline
#else /* !__GNUC__ */
#define INLINE static
#endif /* __GNUC__ */
#endif

/*
 * TETRA run-time optimizations, generally keep these all on, the switches
 * are just there for performance testing
 */

/* define FIXED_SIZE_HEAP to include the fixed size heap code manager, without
   this the allocation are made directly to the resident malloc package;
   see fixheap.[hc] for a description of the fixed size heap manager;
note: this will likely save you LARGE amounts of memory used during
runtime, especially if your system has a power of 2 malloc package
 *this is likely the case*
 */

#if 1 
#define FIXED_SIZE_HEAP
#endif

/* define REG_ARRAY to keep the registers out of the live well; instead
   they will be stored in a separate array that can be indexed much
   quicker than the live well
 */
#ifdef sparc
#undef REG_ARRAY	/* cannot use this with the SPARC implementation */
#else
#define REG_ARRAY	/* cannot use this with the SPARC implementation */
#endif /* sparc */

/* remove this to compare performance to malloc/free allocation */

#include "fixheap.h"

#ifdef FIXED_SIZE_HEAP
    void
ExpandHeap(FixHeapPtr heap)
{
    unsigned long i;

    /* this need to be 'char *'s because they are used for pntr arithmetic */
    char *start, *block;

    if (heap->verbose) {
        fprintf(stdout,
                "## Allocating [%lu] nodes, [%lu] bytes for heap [%s]...\n",
                heap->nodesPerAlloc,
                heap->nodeSize*heap->nodesPerAlloc + sizeof(char *),
                heap->name);
    }

    block =
        (char *)malloc((unsigned int)
                       (heap->nodeSize*heap->nodesPerAlloc + sizeof(char *)));
    if (block == NULL)
    {
        fprintf(stderr, "unable to allocate fixed heap block");
        exit(1);
    }

    start = block;			/* start of allocation */
    block = block + sizeof(char *);	/* start of nodes */

    /* build the free list for this new block */
    for (i=0; i<(heap->nodesPerAlloc-1); i++)
        *((char **)(block+(i*heap->nodeSize))) = 
            (block+((i+1)*heap->nodeSize));

    /* append the new nodes onto the front of the free list */
    *((void **)(block+((heap->nodesPerAlloc-1)*heap->nodeSize))) = 
        heap->freeList;
    heap->freeList = block;

    /* append the block onto the block list */
    *((void **)start) = heap->heapAlloc;
    heap->heapAlloc = start;
}
#endif /* FIXED_SIZE_HEAP */

FixHeapPtr FixHeapCreate(char *name, unsigned long nodeSize, 
        unsigned long nodesPerAlloc, int verbose)
{
    FixHeapPtr heap;

    /* pre-assertions */
    if (nodeSize <= 0 || nodesPerAlloc <= 0)
    {
        fprintf(stderr, "zero size or number");
        abort();
    }

    /* this is required for the free list handling */
    if (nodeSize < sizeof(void *))
    {
        fprintf(stderr, "node size too small");
        abort();
    }

    if (!(heap = (FixHeapPtr)malloc(sizeof(FixHeap))))
    {
        fprintf(stderr, "unable to allocate fixed heap descriptor");
        exit(1);
    }

    heap->name = name;
    heap->nodeSize = nodeSize;
    heap->verbose = verbose;
#ifdef FIXED_SIZE_HEAP
    heap->nodesPerAlloc = nodesPerAlloc;

    /* initially completely empty */
    heap->freeList = NULL;
    heap->heapAlloc = NULL;

    ExpandHeap(heap);	/* initial allocation */
#endif /* FIXED_SIZE_HEAP */

    return heap;
}

void FixHeapRelease(FixHeapPtr heap)
{
    void **block, **nextBlock;

#ifdef FIXED_SIZE_HEAP
    /* free the alloc'ed blocks */
    for (block = (void **)heap->heapAlloc; block != NULL; block = nextBlock) {
        nextBlock = (void **)(*block);
        free(block);
    }
#endif /* FIXED_SIZE_HEAP */

    free(heap);
}

void *FixHeapMalloc(FixHeapPtr heap)
{
    void *node;

#ifdef FIXED_SIZE_HEAP
    if (heap->freeList == NULL)
        ExpandHeap(heap);	/* get more nodes */

    node = heap->freeList;
    heap->freeList = *((void **)heap->freeList);
#else /* !FIXED_SIZE_HEAP */
    node = malloc(heap->nodeSize);
    if (node == NULL)
        printf("unable to allocate heap node");
#endif /* FIXED_SIZE_HEAP */
    return node;
}

void FixHeapFree(FixHeapPtr heap, void *node)
{
#ifdef FIXED_SIZE_HEAP
    *((void **)node) = heap->freeList;
    heap->freeList = node;
#else /* !FIXED_SIZE_HEAP */
    free(node);
#endif /* FIXED_SIZE_HEAP */
}
