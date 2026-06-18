//----------------------------------------------------------------------------
// LAGraph/src/test/test_pref_attach.c: test cases for LAGraph_pref_attach
//----------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

//-----------------------------------------------------------------------------


#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <time.h>

char msg [LAGRAPH_MSG_LEN] ;
void test_pref_attach_1 (void)
{
    GrB_Info info = LAGraph_Init (msg) ;
    TEST_CHECK (info == GrB_SUCCESS) ;

    LAGraph_Graph G = NULL ;
    GrB_Vector Input_i = NULL ;
    GrB_Vector Input_j = NULL ;

    // GrB_Index seed_edges = 10 ;
    // uint64_t seed_i [10] = { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 } ;
    // uint64_t seed_j [10] = { 1, 2, 3, 4, 0, 2, 3, 4, 0, 1 } ;
    GrB_Index seed_edges = 80 ;
    uint64_t seed_i [80] = {
         0,  0,  0,  0,  1,  1,  1,  1,  2,  2,
         2,  2,  3,  3,  3,  3,  4,  4,  4,  4,
         5,  5,  5,  5,  6,  6,  6,  6,  7,  7,
         7,  7,  8,  8,  8,  8,  9,  9,  9,  9,
        10, 10, 10, 10, 11, 11, 11, 11, 12, 12,
        12, 12, 13, 13, 13, 13, 14, 14, 14, 14,
        15, 15, 15, 15, 16, 16, 16, 16, 17, 17,
        17, 17, 18, 18, 18, 18, 19, 19, 19, 19
    } ;
    uint64_t seed_j [80] = {
         1,  5, 10, 15,  2,  6, 11, 16,  3,  7,
        12, 17,  4,  8, 13, 18,  0,  9, 14, 19,
         6, 10, 15,  0,  7, 11, 16,  1,  8, 12,
        17,  2,  9, 13, 18,  3,  0, 14, 19,  4,
        11, 15,  0,  5, 12, 16,  1,  6, 13, 17,
         2,  7, 14, 18,  3,  8, 15, 19,  4,  9,
        16,  0,  5, 10, 17,  1,  6, 11, 18,  2,
         7, 12, 19,  3,  8, 13,  0,  4,  9, 14
    } ;

    OK (GrB_Vector_new (&Input_i, GrB_UINT64, seed_edges)) ;
    OK (GrB_Vector_new (&Input_j, GrB_UINT64, seed_edges)) ;

    uint64_t max_node = 0 ;
    for (GrB_Index k = 0 ; k < seed_edges ; k++)
    {
        OK (GrB_Vector_setElement_UINT64 (Input_i, seed_i [k], k)) ;
        OK (GrB_Vector_setElement_UINT64 (Input_j, seed_j [k], k)) ;
        if (seed_i [k] > max_node) max_node = seed_i [k] ;
        if (seed_j [k] > max_node) max_node = seed_j [k] ;
    }

    uint64_t seed = 42 ;
    GrB_Index batch_size = 5 ;
    double batch_growth_factor = 1.2 ;
    GrB_Index outgoing_edges = 15 ;
    GrB_Index incoming_edges = 15 ;
    bool directed = true ;

    GrB_Index base_node = (GrB_Index) (max_node + 1) ;
    GrB_Index num_nodes = base_node + 1000 ;

    double start = LAGraph_WallClockTime ( ) ;
    info = LAGraph_pref_attach (
        &G,
        num_nodes,
        seed,
        batch_size,
        batch_growth_factor,
        Input_i,
        Input_j,
        incoming_edges,
        outgoing_edges,
        directed,
        msg
    ) ;
    double end = LAGraph_WallClockTime ( ) ;
    double elapsed = end - start ;
    printf("Time taken: %.6f seconds\n", elapsed);

    if (info != GrB_SUCCESS)
    {
        printf ("Error message: %s\n", msg) ;
        GrB_free (&Input_i) ;
        GrB_free (&Input_j) ;
        LAGraph_Finalize (msg) ;
        TEST_CHECK (info == GrB_SUCCESS) ;
        return ;
    }

    TEST_CHECK (G != NULL) ;
    if (G == NULL)
    {
        GrB_free (&Input_i) ;
        GrB_free (&Input_j) ;
        LAGraph_Finalize (msg) ;
        return ;
    }

    OK (LAGraph_Matrix_Print (G->A, LAGraph_SHORT, stdout, msg)) ;

    OK ( LAGraph_Delete (&G, msg)) ;
    OK (GrB_free (&Input_i)) ;
    OK (GrB_free (&Input_j)) ;
    LAGraph_Finalize (msg) ;
}
void test_pref_attach_2 (void)
{
    GrB_Info info = LAGraph_Init (msg) ;
    TEST_CHECK (info == GrB_SUCCESS) ;

    LAGraph_Graph G = NULL ;
    GrB_Vector Input_i = NULL ;
    GrB_Vector Input_j = NULL ;

    // GrB_Index seed_edges = 10 ;
    // uint64_t seed_i [10] = { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 } ;
    // uint64_t seed_j [10] = { 1, 2, 3, 4, 0, 2, 3, 4, 0, 1 } ;
    GrB_Index seed_edges = 80 ;
    uint64_t seed_i [80] = {
         0,  0,  0,  0,  1,  1,  1,  1,  2,  2,
         2,  2,  3,  3,  3,  3,  4,  4,  4,  4,
         5,  5,  5,  5,  6,  6,  6,  6,  7,  7,
         7,  7,  8,  8,  8,  8,  9,  9,  9,  9,
        10, 10, 10, 10, 11, 11, 11, 11, 12, 12,
        12, 12, 13, 13, 13, 13, 14, 14, 14, 14,
        15, 15, 15, 15, 16, 16, 16, 16, 17, 17,
        17, 17, 18, 18, 18, 18, 19, 19, 19, 19
    } ;
    uint64_t seed_j [80] = {
         1,  5, 10, 15,  2,  6, 11, 16,  3,  7,
        12, 17,  4,  8, 13, 18,  0,  9, 14, 19,
         6, 10, 15,  0,  7, 11, 16,  1,  8, 12,
        17,  2,  9, 13, 18,  3,  0, 14, 19,  4,
        11, 15,  0,  5, 12, 16,  1,  6, 13, 17,
         2,  7, 14, 18,  3,  8, 15, 19,  4,  9,
        16,  0,  5, 10, 17,  1,  6, 11, 18,  2,
         7, 12, 19,  3,  8, 13,  0,  4,  9, 14
    } ;

    OK (GrB_Vector_new (&Input_i, GrB_UINT64, seed_edges)) ;
    OK (GrB_Vector_new (&Input_j, GrB_UINT64, seed_edges)) ;

    uint64_t max_node = 0 ;
    for (GrB_Index k = 0 ; k < seed_edges ; k++)
    {
        OK (GrB_Vector_setElement_UINT64 (Input_i, seed_i [k], k)) ;
        OK (GrB_Vector_setElement_UINT64 (Input_j, seed_j [k], k)) ;
        if (seed_i [k] > max_node) max_node = seed_i [k] ;
        if (seed_j [k] > max_node) max_node = seed_j [k] ;
    }

    uint64_t seed = 42 ;
    GrB_Index batch_size = 5 ;
    double batch_growth_factor = 1.2 ;
    GrB_Index outgoing_edges = 15 ;
    GrB_Index incoming_edges = 15 ;
    bool directed = true ;

    GrB_Index base_node = (GrB_Index) (max_node + 1) ;
    GrB_Index num_nodes = base_node + 1000 ;

    info = LAGraph_pref_attach (
        &G,
        num_nodes,
        seed,
        batch_size,
        batch_growth_factor,
        Input_i,
        Input_j,
        incoming_edges,
        outgoing_edges,
        directed,
        msg
    ) ;

    // printf ("\nLAGraph_pref_attach returned: %d\n", info) ;
    if (info != GrB_SUCCESS)
    {
        printf ("Error message: %s\n", msg) ;
        GrB_free (&Input_i) ;
        GrB_free (&Input_j) ;
        LAGraph_Finalize (msg) ;
        TEST_CHECK (info == GrB_SUCCESS) ;
        return ;
    }

    TEST_CHECK (G != NULL) ;
    if (G == NULL)
    {
        GrB_free (&Input_i) ;
        GrB_free (&Input_j) ;
        LAGraph_Finalize (msg) ;
        return ;
    }

    // OK (LAGraph_Matrix_Print (G->A, LAGraph_SHORT, stdout, msg)) ;
    OK (LAGraph_Matrix_Print (G->A, LAGraph_COMPLETE, stdout, msg)) ;

    OK ( LAGraph_Delete (&G, msg)) ;
    OK (GrB_free (&Input_i)) ;
    OK (GrB_free (&Input_j)) ;
    LAGraph_Finalize (msg) ;
}

TEST_LIST =
{
    { "pref_attach_1", test_pref_attach_1 },
    {NULL, NULL}
} ;

