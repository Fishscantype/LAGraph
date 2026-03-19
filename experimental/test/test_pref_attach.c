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

    // int num_nodes,              // total number of nodes in final graph
    // uint64_t seed,              // random seed
    // int batch_size,             // number of new nodes per batch (initial)
    // int batch_growth_factor,    // multiply batch_size by this each iteration
    // GrB_Vector Input_i,         // seed I vector (row indices of seed edges)
    // GrB_Vector Input_j,         // seed J vector (col indices of seed edges)
    // int incoming_edges,         // incoming edges per new node (existing -> new)
    // int outgoing_edges,         // outgoing edges per new node (new -> existing)
    // bool directed,              // (unused for now)
    // double batch_self_edges,    // (unused for now)

    // GrB_Matrix A = NULL ;
    // GrB_Index n = 8 ;
    // GrB_Index nedges = 16 ;
    // double a = 0.45; double b = 0.15; double c = 0.15; double d = 0.25; // RMAT probabilities
    // bool symmetric = false;    // directed graph

    // info = LAGraph_RMAT(&A, n, nedges, a, b, c, d, symmetric, msg);
    // printf("LAGraph_RMAT error: %s\n", msg);
    // TEST_CHECK(info == GrB_SUCCESS);
    // TEST_CHECK(A != NULL);

    // //--------------------------------------------------------------------------
    // // Print the adjacency matrix
    // //--------------------------------------------------------------------------

    // printf("\nAdjacency matrix of RMAT graph:\n");
    // OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

    // //--------------------------------------------------------------------------
    // // free everything and finalize LAGraph
    // //--------------------------------------------------------------------------

    // OK (GrB_free (&A)) ;
    LAGraph_Finalize (msg) ;
}

// void test_RMAT_2 (void)
// {

//     //--------------------------------------------------------------------------
//     // start LAGraph
//     //--------------------------------------------------------------------------

//     GrB_Info info = LAGraph_Init(msg);
//     TEST_CHECK(info == GrB_SUCCESS);

//     //--------------------------------------------------------------------------
//     // Generate a basic RMAT graph
//     //--------------------------------------------------------------------------

//     GrB_Matrix A = NULL ;
//     GrB_Index n = 31 ;           // num of nodes is 2^n
//     GrB_Index nedges = 16 ;     // number of edges
//     double a = 0.45; double b = 0.15; double c = 0.15; double d = 0.25; // RMAT probabilities
//     bool symmetric = false;    // directed graph

//     info = LAGraph_RMAT(&A, n, nedges, a, b, c, d, symmetric, msg);
//     printf("LAGraph_RMAT error: %s\n", msg);
//     TEST_CHECK(info == GrB_SUCCESS);
//     TEST_CHECK(A != NULL);

//     //--------------------------------------------------------------------------
//     // Print the adjacency matrix
//     //--------------------------------------------------------------------------

//     printf("\nAdjacency matrix of RMAT graph:\n");
//     OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

//     //--------------------------------------------------------------------------
//     // free everything and finalize LAGraph
//     //--------------------------------------------------------------------------

//     OK (GrB_free (&A)) ;
//     LAGraph_Finalize (msg) ;
// }

// void test_RMAT_3 (void)
// {

//     //--------------------------------------------------------------------------
//     // start LAGraph
//     //--------------------------------------------------------------------------

//     GrB_Info info = LAGraph_Init(msg);
//     TEST_CHECK(info == GrB_SUCCESS);

//     //--------------------------------------------------------------------------
//     // Generate a basic RMAT graph
//     //--------------------------------------------------------------------------

//     GrB_Matrix A = NULL ;
//     GrB_Index n = 32 ;           // num of nodes is 2^n
//     GrB_Index nedges = 16 ;     // number of edges
//     double a = 0.45; double b = 0.15; double c = 0.15; double d = 0.25; // RMAT probabilities
//     bool symmetric = false;    // directed graph

//     info = LAGraph_RMAT(&A, n, nedges, a, b, c, d, symmetric, msg);
//     printf("LAGraph_RMAT error: %s\n", msg);
//     TEST_CHECK(info == GrB_SUCCESS);
//     TEST_CHECK(A != NULL);

//     //--------------------------------------------------------------------------
//     // Print the adjacency matrix
//     //--------------------------------------------------------------------------

//     printf("\nAdjacency matrix of RMAT graph:\n");
//     OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

//     //--------------------------------------------------------------------------
//     // free everything and finalize LAGraph
//     //--------------------------------------------------------------------------

//     OK (GrB_free (&A)) ;
//     LAGraph_Finalize (msg) ;
// }

// void test_RMAT_4 (void)
// {

//     //--------------------------------------------------------------------------
//     // start LAGraph
//     //--------------------------------------------------------------------------

//     GrB_Info info = LAGraph_Init(msg);
//     TEST_CHECK(info == GrB_SUCCESS);

//     //--------------------------------------------------------------------------
//     // Generate a basic RMAT graph
//     //--------------------------------------------------------------------------

//     GrB_Matrix A = NULL ;
//     GrB_Index n = 61 ;           // num of nodes is 2^n
//     GrB_Index nedges = 16 ;     // number of edges
//     double a = 0.45; double b = 0.15; double c = 0.15; double d = 0.25; // RMAT probabilities
//     bool symmetric = false;    // directed graph

//     info = LAGraph_RMAT(&A, n, nedges, a, b, c, d, symmetric, msg);
//     printf("LAGraph_RMAT error: %s\n", msg);
//     TEST_CHECK(info == GrB_SUCCESS);
//     TEST_CHECK(A != NULL);

//     //--------------------------------------------------------------------------
//     // Print the adjacency matrix
//     //--------------------------------------------------------------------------

//     printf("\nAdjacency matrix of RMAT graph:\n");
//     OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

//     //--------------------------------------------------------------------------
//     // free everything and finalize LAGraph
//     //--------------------------------------------------------------------------

//     OK (GrB_free (&A)) ;
//     LAGraph_Finalize (msg) ;
// }

// void test_RMAT_5 (void)
// {
//     // Large graph, timed

//     GrB_Info info = LAGraph_Init(msg);
//     TEST_CHECK(info == GrB_SUCCESS);


//     GrB_Matrix A = NULL ;
//     GrB_Index n = 59 ;           // num of nodes is 2^n
//     GrB_Index nedges = 100000000 ;     // number of edges
//     double a = 0.45; double b = 0.15; double c = 0.15; double d = 0.25; // RMAT probabilities
//     bool symmetric = false;    // directed graph

//     // clock_t start, end;
//     // double cpu_time_used;
//     // start = clock();
//     struct timespec start, end;
//     clock_gettime(CLOCK_MONOTONIC, &start);
//     info = LAGraph_RMAT(&A, n, nedges, a, b, c, d, symmetric, msg);
//     clock_gettime(CLOCK_MONOTONIC, &end);
//     double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
//     printf("Time taken: %.6f seconds\n", elapsed);
//     // printf("LAGraph_RMAT error: %s\n", msg);
//     TEST_CHECK(info == GrB_SUCCESS);
//     TEST_CHECK(A != NULL);

//     //--------------------------------------------------------------------------
//     // Print the adjacency matrix
//     //--------------------------------------------------------------------------

//     // printf("\nAdjacency matrix of RMAT graph:\n");
//     // OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

//     // end = clock();
//     // cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
//     // printf("Time taken: %f seconds\n", cpu_time_used);


//     //--------------------------------------------------------------------------
//     // free everything and finalize LAGraph
//     //--------------------------------------------------------------------------

//     OK (GrB_free (&A)) ;
//     LAGraph_Finalize (msg) ;
// }

// void test_RMAT_6 (void)
// {
//     // Large graph, print output to file

//     //--------------------------------------------------------------------------
//     // start LAGraph
//     //--------------------------------------------------------------------------

//     GrB_Info info = LAGraph_Init(msg);
//     // TEST_CHECK(info == GrB_SUCCESS);

//     //--------------------------------------------------------------------------
//     // Generate a basic RMAT graph
//     //--------------------------------------------------------------------------

//     GrB_Matrix A = NULL ;
//     GrB_Index n = 16 ;           // num of nodes is 2^n
//     GrB_Index nedges = 100000 ;     // number of edges
//     double a = 0.25; double b = 0.25; double c = 0.25; double d = 0.25; // RMAT probabilities
//     bool symmetric = false;    // d`irected graph

//     info = LAGraph_RMAT(&A, n, nedges, a, b, c, d, symmetric, msg);
//     // TEST_CHECK(info == GrB_SUCCESS);
//     // TEST_CHECK(A != NULL);

//     //--------------------------------------------------------------------------
//     // Print the adjacency matrix
//     //--------------------------------------------------------------------------

//     OK(LAGraph_Matrix_Print(A, LAGraph_COMPLETE, stdout, msg));

//     //--------------------------------------------------------------------------
//     // free everything and finalize LAGraph
//     //--------------------------------------------------------------------------

//     OK (GrB_free (&A)) ;
//     LAGraph_Finalize (msg) ;
// }

TEST_LIST =
{
    { "pref_attach_1", test_pref_attach_1 },
    // {"RMAT_Generate_basic", test_RMAT_1},
    // {"RMAT_Generate_31", test_RMAT_2},
    // {"RMAT_Generate_32", test_RMAT_3},
    // {"RMAT_Generate_62", test_RMAT_4},
    // {"RMAT_Generate_large", test_RMAT_5},
    // {"RMAT_Generate_large_print", test_RMAT_6},
    {NULL, NULL}
} ;

