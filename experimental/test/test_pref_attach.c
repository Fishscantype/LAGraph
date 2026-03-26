//----------------------------------------------------------------------------
// LAGraph/src/test/test_HelloWorld.c: test cases for LAGraph_HelloWorld
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
    GrB_Info info = LAGraph_Init(msg);
    TEST_CHECK(info == GrB_SUCCESS);

    GrB_Matrix A = NULL ;
    GrB_Vector Input_i = NULL ;
    GrB_Vector Input_j = NULL ;

    GrB_Index seed_edges = 3 ;
    uint64_t seed_i[] = { 0, 1, 2 } ;
    uint64_t seed_j[] = { 1, 2, 0 } ;

    OK (GrB_Vector_new (&Input_i, GrB_UINT64, seed_edges)) ;
    OK (GrB_Vector_new (&Input_j, GrB_UINT64, seed_edges)) ;

    for (GrB_Index k = 0 ; k < seed_edges ; k++)
    {
        OK (GrB_Vector_setElement_UINT64 (Input_i, seed_i[k], k)) ;
        OK (GrB_Vector_setElement_UINT64 (Input_j, seed_j[k], k)) ;
    }

    printf ("\n--- Seed Input_i ---\n") ;
    OK (LAGraph_Vector_Print (Input_i, LAGraph_COMPLETE, stdout, msg)) ;
    printf ("\n--- Seed Input_j ---\n") ;
    OK (LAGraph_Vector_Print (Input_j, LAGraph_COMPLETE, stdout, msg)) ;

    GrB_Index num_nodes = 10 ;
    uint64_t seed = 42 ;
    GrB_Index batch_size = 5 ;
    double batch_growth_factor = 1.0 ;
    GrB_Index outgoing_edges = 2 ;
    GrB_Index incoming_edges = 0 ;
    bool directed = false ;
    double batch_alpha = 0.0 ;

    info = LAGraph_pref_attach (
        &A,
        num_nodes,
        seed,
        batch_size,
        batch_growth_factor,
        Input_i,
        Input_j,
        incoming_edges,
        outgoing_edges,
        directed,
        batch_alpha,
        msg
    ) ;

    printf ("\nLAGraph_pref_attach returned: %d\n", info) ;
    if (info != GrB_SUCCESS){
        printf ("Error message: %s\n", msg) ;
        GrB_free (&Input_i) ;
        GrB_free (&Input_j) ;
        LAGraph_Finalize (msg) ;
        TEST_CHECK (info == GrB_SUCCESS) ;
        return ;
    }
    TEST_CHECK (A != NULL) ;

    printf ("\n--- Output adjacency matrix ---\n") ;
    OK (LAGraph_Matrix_Print (A, LAGraph_COMPLETE, stdout, msg)) ;

    GrB_Index nrows, ncols, nvals ;
    OK (GrB_Matrix_nrows (&nrows, A)) ;
    OK (GrB_Matrix_ncols (&ncols, A)) ;
    OK (GrB_Matrix_nvals (&nvals, A)) ;
    printf ("\nMatrix: %lu x %lu, nvals = %lu\n",
        (unsigned long) nrows, (unsigned long) ncols,
        (unsigned long) nvals) ;

    OK (GrB_free (&A)) ;
    OK (GrB_free (&Input_i)) ;
    OK (GrB_free (&Input_j)) ;
    LAGraph_Finalize (msg) ;
}

void test_pref_attach_2 (void)
{
    GrB_Info info = LAGraph_Init (msg) ;
    TEST_CHECK (info == GrB_SUCCESS) ;

    GrB_Matrix A = NULL ;
    GrB_Vector Input_i = NULL ;
    GrB_Vector Input_j = NULL ;

    GrB_Index seed_edges = 10 ;
    uint64_t seed_i [10] = { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 } ;
    uint64_t seed_j [10] = { 1, 2, 3, 4, 0, 2, 3, 4, 0, 1 } ;

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
    GrB_Index incoming_edges = 1 ;
    bool directed = true ;
    double batch_alpha = 1.0 ;

    GrB_Index base_node = (GrB_Index) (max_node + 1) ;
    GrB_Index num_nodes = base_node + 1000 ;

    info = LAGraph_pref_attach (
        &A,
        num_nodes,
        seed,
        batch_size,
        batch_growth_factor,
        Input_i,
        Input_j,
        incoming_edges,
        outgoing_edges,
        directed,
        batch_alpha,
        msg
    ) ;

    TEST_CHECK (info == GrB_SUCCESS) ;
    TEST_CHECK (A != NULL) ;

    OK (LAGraph_Matrix_Print (A, LAGraph_SHORT, stdout, msg)) ;
    // OK (LAGraph_Matrix_Print (A, LAGraph_COMPLETE, stdout, msg)) ;

    OK (GrB_free (&A)) ;
    OK (GrB_free (&Input_i)) ;
    OK (GrB_free (&Input_j)) ;
    LAGraph_Finalize (msg) ;
}

TEST_LIST =
{
    // { "pref_attach_1", test_pref_attach_1 },
    { "pref_attach_2", test_pref_attach_2 },
    {NULL, NULL}
} ;

