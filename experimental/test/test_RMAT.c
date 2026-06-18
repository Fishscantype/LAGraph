//----------------------------------------------------------------------------
// LAGraph/src/test/test_RMAT.c: test cases for LAGraph_RMAT
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

static void run_rmat_test
(
    GrB_Index n,
    GrB_Index nedges,
    double a,
    double b,
    double c,
    double d,
    uint64_t seed,
    bool print_graph
)
{
    GrB_Info info = LAGraph_Init (msg) ;
    TEST_CHECK (info == GrB_SUCCESS) ;

    LAGraph_Graph G = NULL ;

    info = LAGraph_RMAT (&G, n, nedges, seed, a, b, c, d, msg) ;
    TEST_CHECK (info == GrB_SUCCESS) ;
    TEST_CHECK (G != NULL) ;

    if (print_graph)
    {
        printf ("\nRMAT graph:\n") ;
        OK (LAGraph_Graph_Print (G, LAGraph_COMPLETE, stdout, msg)) ;
    }

    OK (LAGraph_Delete (&G, msg)) ;
    LAGraph_Finalize (msg) ;
}

void test_RMAT_1 (void)
{
    run_rmat_test (8, 16, 0.45, 0.15, 0.15, 0.25, 42, true) ;
}

void test_RMAT_2 (void)
{
    run_rmat_test (31, 16, 0.45, 0.15, 0.15, 0.25, 42, true) ;
}

void test_RMAT_3 (void)
{
    run_rmat_test (32, 16, 0.45, 0.15, 0.15, 0.25, 42, true) ;
}

void test_RMAT_4 (void)
{
    run_rmat_test (60, 16, 0.45, 0.15, 0.15, 0.25, 42, true) ;
}

void test_RMAT_5 (void)
{
    // Large graph, timed

    GrB_Info info = LAGraph_Init (msg) ;
    TEST_CHECK (info == GrB_SUCCESS) ;

    LAGraph_Graph G = NULL ;

    double start = LAGraph_WallClockTime ( ) ;
    info = LAGraph_RMAT (&G, 20, 100000, 42, 0.5, 0.1, 0.1, 0.3, msg) ;
    double end = LAGraph_WallClockTime ( ) ;
    double elapsed = end - start ;
    printf ("Time taken: %.6f seconds\n", elapsed) ;

    TEST_CHECK (info == GrB_SUCCESS) ;
    TEST_CHECK (G != NULL) ;

    printf ("\nRMAT graph:\n") ;
    OK (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;

    OK (LAGraph_Delete (&G, msg)) ;
    LAGraph_Finalize (msg) ;
}

TEST_LIST =
{
    {"RMAT_Generate_basic", test_RMAT_1},
    {"RMAT_Generate_31", test_RMAT_2},
    {"RMAT_Generate_32", test_RMAT_3},
    {"RMAT_Generate_60", test_RMAT_4},
    {"RMAT_Generate_large", test_RMAT_5},
    {NULL, NULL}
} ;

