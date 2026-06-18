//------------------------------------------------------------------------------
// LAGraph_Random_Graph: generate a random graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Matthew Fisher, Texas A&M University

//------------------------------------------------------------------------------

#define LG_FREE_WORK  \
    {                 \
        GrB_free(&A); \
    }

#define LG_FREE_ALL             \
    {                           \
        LG_FREE_WORK;           \
        LAGraph_Delete(G, msg); \
    }

#include "LG_internal.h"
#include "LAGraphX.h"

//------------------------------------------------------------------------------
// LAGraph_Random_Graph
//------------------------------------------------------------------------------

GrB_Info LAGraph_Random_Graph // random graph of any built-in type
    (
        // output
        LAGraph_Graph *G,      // G is constructed on output
        // input
        GrB_Type type,         // type of matrix to construct
        LAGraph_Kind kind,     // type of graph
        bool allow_self_edges, // should this graph allow self edges
        bool weighted,         // should this graph be weighted
        GrB_Index nnodes,      // # of nodes in G
        int64_t nedges,        // approximate # of edges in G
        uint64_t seed,         // random number seed
        char *msg)
{
    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    GrB_Matrix A = NULL;

    LG_ASSERT(G != NULL, GrB_NULL_POINTER);
    LG_ASSERT(type != NULL, GrB_NULL_POINTER);
    LG_ASSERT(kind != LAGraph_KIND_UNKNOWN, GrB_INVALID_VALUE);
    LG_ASSERT(nnodes > 0, GrB_INVALID_VALUE);
    LG_ASSERT(nedges > 0, GrB_INVALID_VALUE);


    //--------------------------------------------------------------------------
    // Initialize values
    //--------------------------------------------------------------------------

    (*G) = NULL;

    double density = (double) nedges / (nnodes * nnodes);

    GRB_TRY(LAGraph_Random_Matrix(&A, type, nnodes, nnodes, density, seed, msg));

    //--------------------------------------------------------------------------
    // Make symmetric if undirected graph is desired
    //--------------------------------------------------------------------------

    if (kind == LAGraph_ADJACENCY_UNDIRECTED)
    {
        GrB_Matrix L = NULL, LT = NULL, Out = NULL;

        GRB_TRY(GrB_Matrix_new(&L, type, nnodes, nnodes));
        GRB_TRY(GrB_select(L, NULL, NULL, GrB_TRIL, A, (int64_t)(-1), NULL));

        GRB_TRY(GrB_Matrix_new(&LT, type, nnodes, nnodes));
        GRB_TRY(GrB_transpose(LT, NULL, NULL, L, NULL));

        GRB_TRY(GrB_Matrix_new(&Out, type, nnodes, nnodes));
        GrB_BinaryOp operator = weighted ? GrB_PLUS_FP64 : GrB_ONEB_FP64;
        GRB_TRY(GrB_eWiseAdd(Out, NULL, NULL, operator, L, LT, NULL));

        GrB_free(&A);
        GrB_free(&L);
        GrB_free(&LT);
        A = Out;
    }

    if (!weighted)
    {
        GRB_TRY(GrB_assign(A, A, NULL, 1, GrB_ALL, nnodes, GrB_ALL, nnodes, GrB_DESC_S));
    }

    //--------------------------------------------------------------------------
    // Create graph from matrix
    //--------------------------------------------------------------------------

    LAGRAPH_TRY(LAGraph_New(G, &A, kind, msg));
    if (!allow_self_edges)
    {
        LAGRAPH_TRY(LAGraph_DeleteSelfEdges(*G, msg));
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LG_FREE_WORK;
    return (GrB_SUCCESS);
}
