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

// This is a very simple "hello world" example of a test program for an
// algorithm in the experimental/algorithm folder.

// TODO: GrB_select with GrB_TRIL to grab the lower triangular matric. Do -1 to get rid of the diagonal to prevent self edges
// make random matrix wrapper to create random graph, make graph square, then use above step to make symmetric IF directed is false
// operator for ewise_add does not matter if tril and tril transpose are added, can just pass in grb_fp64 as type
// operator: GrB_oneb_fp64 for unweighted graph
// function for removing self edges: LAGrapg_delete_self_edges

#include <stdio.h>
#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>
#include <time.h>

char msg [LAGRAPH_MSG_LEN] ;

void test_RMAT_1 (void)
{

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    GrB_Info info = LAGraph_Init(msg);
    TEST_CHECK(info == GrB_SUCCESS);

    //--------------------------------------------------------------------------
    // Generate a basic RMAT graph
    //--------------------------------------------------------------------------

    GrB_Matrix A = NULL ;
    GrB_Index n = 8 ;           // num of nodes is 2^n
    GrB_Index nedges = 16 ;     // number of edges
    double a = 0.45; double b = 0.15; double c = 0.15; double d = 0.25; // RMAT probabilities
    uint64_t seed = 42 ;

    info = LAGraph_RMAT(&A, n, nedges, seed, a, b, c, d, msg);
    TEST_CHECK(info == GrB_SUCCESS);
    TEST_CHECK(A != NULL);

    //--------------------------------------------------------------------------
    // Print the adjacency matrix
    //--------------------------------------------------------------------------

    printf("\nAdjacency matrix of RMAT graph:\n");
    OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------

    OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
}

void test_RMAT_2 (void)
{

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    GrB_Info info = LAGraph_Init(msg);
    TEST_CHECK(info == GrB_SUCCESS);

    //--------------------------------------------------------------------------
    // Generate a basic RMAT graph
    //--------------------------------------------------------------------------

    GrB_Matrix A = NULL ;
    GrB_Index n = 31 ;           // num of nodes is 2^n
    GrB_Index nedges = 16 ;     // number of edges
    double a = 0.45; double b = 0.15; double c = 0.15; double d = 0.25; // RMAT probabilities
    uint64_t seed = 42 ;

    info = LAGraph_RMAT(&A, n, nedges, seed, a, b, c, d, msg);
    TEST_CHECK(info == GrB_SUCCESS);
    TEST_CHECK(A != NULL);

    //--------------------------------------------------------------------------
    // Print the adjacency matrix
    //--------------------------------------------------------------------------

    printf("\nAdjacency matrix of RMAT graph:\n");
    OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------

    OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
}

void test_RMAT_3 (void)
{

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    GrB_Info info = LAGraph_Init(msg);
    TEST_CHECK(info == GrB_SUCCESS);

    //--------------------------------------------------------------------------
    // Generate a basic RMAT graph
    //--------------------------------------------------------------------------

    GrB_Matrix A = NULL ;
    GrB_Index n = 32 ;           // num of nodes is 2^n
    GrB_Index nedges = 16 ;     // number of edges
    double a = 0.45; double b = 0.15; double c = 0.15; double d = 0.25; // RMAT probabilities
    uint64_t seed = 42 ;

    info = LAGraph_RMAT(&A, n, nedges, seed, a, b, c, d, msg);
    TEST_CHECK(info == GrB_SUCCESS);
    TEST_CHECK(A != NULL);

    //--------------------------------------------------------------------------
    // Print the adjacency matrix
    //--------------------------------------------------------------------------

    printf("\nAdjacency matrix of RMAT graph:\n");
    OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------

    OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
}

void test_RMAT_4 (void)
{

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    GrB_Info info = LAGraph_Init(msg);
    TEST_CHECK(info == GrB_SUCCESS);

    //--------------------------------------------------------------------------
    // Generate a basic RMAT graph
    //--------------------------------------------------------------------------

    GrB_Matrix A = NULL ;
    GrB_Index n = 60 ;           // num of nodes is 2^n
    GrB_Index nedges = 16 ;     // number of edges
    double a = 0.45; double b = 0.15; double c = 0.15; double d = 0.25; // RMAT probabilities
    uint64_t seed = 42 ;

    info = LAGraph_RMAT(&A, n, nedges, seed, a, b, c, d, msg);
    TEST_CHECK(info == GrB_SUCCESS);
    TEST_CHECK(A != NULL);

    //--------------------------------------------------------------------------
    // Print the adjacency matrix
    //--------------------------------------------------------------------------

    printf("\nAdjacency matrix of RMAT graph:\n");
    OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------

    OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
}

void test_RMAT_5 (void)
{
    // Large graph, timed

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    GrB_Info info = LAGraph_Init(msg);
    TEST_CHECK(info == GrB_SUCCESS);

    //--------------------------------------------------------------------------
    // Generate a basic RMAT graph
    //--------------------------------------------------------------------------

    GrB_Matrix A = NULL ;
    GrB_Index n = 59 ;           // num of nodes is 2^n
    GrB_Index nedges = 100000000 ;     // number of edges
    double a = 0.45; double b = 0.15; double c = 0.15; double d = 0.25; // RMAT probabilities
    uint64_t seed = 42 ;

    // clock_t start, end;
    // double cpu_time_used;
    // start = clock();
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    info = LAGraph_RMAT(&A, n, nedges, seed, a, b, c, d, msg);
    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Time taken: %.6f seconds\n", elapsed);
    // printf("LAGraph_RMAT error: %s\n", msg);
    TEST_CHECK(info == GrB_SUCCESS);
    TEST_CHECK(A != NULL);

    //--------------------------------------------------------------------------
    // Print the adjacency matrix
    //--------------------------------------------------------------------------

    // printf("\nAdjacency matrix of RMAT graph:\n");
    // OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

    // end = clock();
    // cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    // printf("Time taken: %f seconds\n", cpu_time_used);


    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------

    OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
}

void test_RMAT_6 (void)
{
    // Large graph, print output to file

    //--------------------------------------------------------------------------
    // start LAGraph
    //--------------------------------------------------------------------------

    GrB_Info info = LAGraph_Init(msg);
    // TEST_CHECK(info == GrB_SUCCESS);

    //--------------------------------------------------------------------------
    // Generate a basic RMAT graph
    //--------------------------------------------------------------------------

    GrB_Matrix A = NULL ;
    GrB_Index n = 16 ;           // num of nodes is 2^n
    GrB_Index nedges = 100000 ;     // number of edges
    double a = 0.25; double b = 0.25; double c = 0.25; double d = 0.25; // RMAT probabilities
    uint64_t seed = 42 ;

    info = LAGraph_RMAT(&A, n, nedges, seed, a, b, c, d, msg);
    // TEST_CHECK(info == GrB_SUCCESS);
    // TEST_CHECK(A != NULL);

    //--------------------------------------------------------------------------
    // Print the adjacency matrix
    //--------------------------------------------------------------------------

    OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

    //--------------------------------------------------------------------------
    // free everything and finalize LAGraph
    //--------------------------------------------------------------------------

    OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
}

TEST_LIST =
{
    {"RMAT_Generate_basic", test_RMAT_1},
    {"RMAT_Generate_31", test_RMAT_2},
    {"RMAT_Generate_32", test_RMAT_3},
    {"RMAT_Generate_62", test_RMAT_4},
    {"RMAT_Generate_large", test_RMAT_5},
    // {"RMAT_Generate_large_print", test_RMAT_6},
    {NULL, NULL}
} ;

