//------------------------------------------------------------------------------
// LAGraph_pref_attach: an implementation of the preferential attachment graph generation prpcess
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

// TODO: intro section

#define LG_FREE_WORK                            \
{                                               \
    GrB_free (&State) ;                         \
    GrB_Vector_free (&Scaled_state) ;           \
    GrB_Vector_free (&Output_i) ;               \
    GrB_Vector_free (&Output_j) ;               \
    GrB_Scalar_free (&Scalar_one) ;             \
    GrB_Vector_free (&Batch_random_in) ;        \
    GrB_Vector_free (&Batch_random_out) ;       \
    GrB_Vector_free (&Batch_gather) ;           \
    GrB_Vector_free (&Batch_chunk) ;            \
    GrB_IndexUnaryOp_free (&scale_op) ;         \
    GrB_IndexUnaryOp_free (&chunk_index_op) ;   \
    GrB_Type_free (&scale_param_type) ;         \
    GrB_Type_free (&chunk_index_param_type) ;   \
    GrB_Scalar_free (&scale_scalar) ;           \
    GrB_Scalar_free (&chunk_index_scalar) ;     \
}

#define LG_FREE_ALL                         \
{                                           \
    LG_FREE_WORK ;                          \
    GrB_free (&Y) ;                         \
}

#include "LG_internal.h"
#include "LAGraphX.h"

typedef struct
{
    uint64_t base_edge ;
} Scale_Params ;
#define SCALE_PARAMS                \
"typedef struct {                   \n"\
"    uint64_t base_edge ;           \n"\
"} Scale_Params ;                   "

typedef struct
{
    uint64_t base_node ;
    uint64_t edges ;
} Chunk_Index_Params ;

#define CHUNK_INDEX_PARAMS          \
"typedef struct {                   \n"\
"    uint64_t base_node ;           \n"\
"    uint64_t edges ;               \n"\
"} Chunk_Index_Params ;             "

void pref_attach_scale_op
(
    void *z,
    const void *x,
    GrB_Index i,
    GrB_Index j,
    const void *y
)
{
    const uint64_t *r = (const uint64_t *) x ;
    const Scale_Params *p = (const Scale_Params *) y ;
    if (p->base_edge == 0){
        *((uint64_t *) z) = 0 ;
    } else {
        double scaled = ((double)(*r) / (double) UINT64_MAX)
            * (double)(p->base_edge - 1) ;
        *((uint64_t *) z) = (uint64_t) scaled ;
    }
}
#define PREF_ATTACH_SCALE_OP                                              \
"void pref_attach_scale_op(                                            \n"\
"    void *z,                                                            \n"\
"    const void *x,                                                      \n"\
"    GrB_Index i,                                                        \n"\
"    GrB_Index j,                                                 \n"\
"    const void *y                                                       \n"\
")                                                                       \n"\
"{                                                                       \n"\
"    const uint64_t *r = (const uint64_t *) x ;                   \n"\
"    const Scale_Params *p = (const Scale_Params *) y ;                  \n"\
"    if (p->base_edge == 0){                                            \n"\
"        *((uint64_t *) z) = 0 ;                                        \n"\
"    } else {                                                            \n"\
"        double scaled = ((double)(*r) / (double) UINT64_MAX)    \n"\
"            * (double)(p->base_edge - 1)) ;                              \n"\
"        *((uint64_t *) z) = (uint64_t) scaled ;                        \n"\
"    }                                                                   \n"\
"}                                                                       "

void pref_attach_chunk_index_op
(
    void *z,
    const void *x,
    GrB_Index i,
    GrB_Index j,
    const void *y
)
{
    const Chunk_Index_Params *p = (const Chunk_Index_Params *) y ;
    *((uint64_t *) z) = p->base_node + (uint64_t) i / p->edges ;
}

#define PREF_ATTACH_CHUNK_INDEX_OP                                               \
"void pref_attach_chunk_index_op(                                             \n"\
"    void *z,                                                            \n"\
"    const void *x,                                                      \n"\
"    GrB_Index i,                                                        \n"\
"    GrB_Index j,                                                 \n"\
"    const void *y                                                       \n"\
")                                                                       \n"\
"{                                                                       \n"\
"    const Chunk_Index_Params *p = (const Chunk_Index_Params *) y ;                    \n"\
"    *((uint64_t *) z) = p->base_node + (uint64_t) i / p->edges;            \n"\
"}                                                                       "

int LAGraph_pref_attach
(
    // output
    GrB_Matrix *Yhandle,    // Y, created on output
    // input:
    int num_nodes,
    uint64_t seed,
    int batch_size,
    int batch_growth_factor,
    GrB_Vector Input_i,
    GrB_Vector Input_j,
    int incoming_edges,
    int outgoing_edges,
    bool directed, // if directed, do outgoing edges and mirror
    double batch_self_edges, // mult. factor for batch corner density calculation
    char *msg
)
{

    //--------------------------------------------------------------------------
    // Set up workspace
    //--------------------------------------------------------------------------

    GrB_Vector State = NULL ;
    GrB_Vector Scaled_state = NULL ;
    GrB_Vector Output_i = NULL ;
    GrB_Vector Output_j = NULL ;
    GrB_Scalar Scalar_one = NULL ;
    GrB_Vector Batch_random_in = NULL ;
    GrB_Vector Batch_random_out = NULL ;
    GrB_Vector Batch_gather = NULL ;
    GrB_Vector Batch_chunk = NULL ;
    GrB_IndexUnaryOp scale_op = NULL ;
    GrB_IndexUnaryOp chunk_index_op = NULL ;
    GrB_Type scale_param_type = NULL ;
    GrB_Type chunk_index_param_type = NULL ;
    GrB_Scalar scale_scalar = NULL ;
    GrB_Scalar chunk_index_scalar = NULL ;
    GrB_Matrix Y = NULL ;
    
    LG_CLEAR_MSG ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_ASSERT (Yhandle != NULL, GrB_NULL_POINTER) ;
    (*Yhandle) = NULL ;

    LG_ASSERT (Input_i != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (Input_j != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT_MSG (num_nodes > 0, GrB_INVALID_VALUE, "num_nodes must be positive") ;
    LG_ASSERT_MSG (outgoing_edges > 0, GrB_INVALID_VALUE, "outgoing_edges must be positive") ;
    LG_ASSERT_MSG (incoming_edges >= 0, GrB_INVALID_VALUE, "incoming_edges must be non-negative") ;
    LG_ASSERT_MSG (batch_size > 0, GrB_INVALID_VALUE, "batch_size must be positive") ;
    LG_ASSERT_MSG (batch_growth_factor >= 1, GrB_INVALID_VALUE, "batch_growth_factor must be >= 1") ;
 
    // //--------------------------------------------------------------------------
    // // Evaluate seed vectors
    // //--------------------------------------------------------------------------

    // GrB_Index seed_nvals_i ;
    // GRB_TRY (GrB_Vector_nvals (&seed_nvals_i, Input_i)) ;
    // GrB_Index seed_nvals_j ;
    // GRB_TRY (GrB_Vector_nvals (&seed_nvals_j, Input_j)) ;
    // LG_ASSERT_MSG (seed_nvals_i == seed_nvals_j, GrB_DIMENSION_MISMATCH, "Input vectors must have same number of edges") ;

    // //TODO: get max node value from I and J
    // uint64_t max_node = seed_nvals_i * 2 ;
    // uint64_t base_node = max_node + 1 ;

    // // //--------------------------------------------------------------------------
    // // // Calculate num_edges
    // // //--------------------------------------------------------------------------

    // if(!directed){
    //     incoming_edges = 0 ;
    // }
    // uint64_t num_nodes_to_add = (num_nodes < base_node ? 0 : base_node - num_nodes) ;
    // uint64_t num_edges = num_nodes_to_add * (incoming_edges + outgoing_edges) + seed_nvals_i ;


    // // //--------------------------------------------------------------------------
    // // // Calculate num_batches
    // // //--------------------------------------------------------------------------

    // uint64_t num_batches = 0 ;
    // uint64_t nodes_counted = 0 ;
    // uint64_t b_size = batch_size ;
    // while(nodes_counted < num_nodes_to_add){
    //     uint64_t temp_size = b_size ;
    //     if (nodes_counted + temp_size > num_nodes_to_add){
    //         temp_size = num_nodes_to_add - nodes_counted ;
    //     }
    //     nodes_counted += temp_size ;
    //     ++num_batches ;
    //     b_size *= batch_growth_factor ;
    // }
    
    // // //--------------------------------------------------------------------------
    // // // Initialize Output_i and Output_j vectors
    // // //--------------------------------------------------------------------------

    // //TODO: come back and reduce to uint32_t if possible
    // GRB_TRY (GrB_Vector_new (&Output_i, GrB_UINT64, num_edges)) ;
    // GRB_TRY (GrB_Vector_new (&Output_j, GrB_UINT64, num_edges)) ;

    // GrB_Index range_seed[3] ;
    // range_seed[GxB_BEGIN] = 0 ;
    // range_seed[GxB_END] = seed_nvals_i - 1 ;
    // range_seed[GxB_INC] = 1 ;

    // GRB_TRY (GxB_Vector_assign (Output_i, NULL, NULL, Input_i, range_seed, GxB_RANGE, NULL)) ;
    // GRB_TRY (GxB_Vector_assign (Output_j, NULL, NULL, Input_j, range_seed, GxB_RANGE, NULL)) ;

    // // //--------------------------------------------------------------------------
    // // // Create random state vector
    // // //--------------------------------------------------------------------------

    // GRB_TRY (GrB_Vector_new (&State, GrB_UINT64, num_edges)) ;
    // GRB_TRY (GrB_assign (State, NULL, NULL, 0, GrB_ALL, num_edges, NULL)) ;
    // LG_TRY (LAGraph_Random_Seed (State, seed, msg)) ;

    // // //--------------------------------------------------------------------------
    // // // Prepare parameter structs, operators, and workspace vectors
    // // //--------------------------------------------------------------------------

    // GRB_TRY (GrB_Type_new (&scale_param_type, sizeof (Scale_Params))) ;
    // GRB_TRY (GrB_Type_new (&chunk_index_param_type, sizeof (Chunk_Index_Params))) ;

    // GRB_TRY (GrB_Scalar_new (&scale_scalar, scale_param_type)) ;
    // GRB_TRY (GrB_Scalar_new (&chunk_index_scalar, chunk_index_param_type)) ;

    // //TODO: come back and reduce to uint32_t if possible
    // GRB_TRY (GrB_IndexUnaryOp_new (&scale_op, (GxB_index_unary_function) pref_attach_scale_op, GrB_UINT64, GrB_UINT64, scale_param_type)) ;
    // GRB_TRY (GrB_IndexUnaryOp_new (&chunk_index_op, (GxB_index_unary_function) pref_attach_chunk_index_op, GrB_UINT64, GrB_UINT64, chunk_index_param_type)) ;

    // GRB_TRY (GrB_Vector_new (&Scaled_state, GrB_UINT64, num_edges)) ;
    // GRB_TRY (GrB_Vector_new (&Batch_random_out, GrB_UINT64, num_edges)) ;
    // GRB_TRY (GrB_Vector_new (&Batch_random_in, GrB_UINT64, num_edges)) ;
    // GRB_TRY (GrB_Vector_new (&Batch_gather, GrB_UINT64, num_edges)) ;
    // GRB_TRY (GrB_Vector_new (&Batch_chunk, GrB_UINT64, num_edges)) ;

    // // //--------------------------------------------------------------------------
    // // // Begin main loop
    // // //--------------------------------------------------------------------------

    // uint64_t curr_node = base_node ;
    // GrB_Index curr_edge = seed_nvals_i ;
    // uint64_t curr_batch_size = batch_size ;

    // for(int batch = 0; batch < num_batches; ++batch){
    //     // Prepare variables for batch
    //     uint64_t nodes_in_batch = (curr_batch_size - 1 > num_nodes - curr_node ? num_nodes - curr_node : curr_batch_size) ;
    //     GrB_Index incoming_edges_in_batch = (GrB_Index) nodes_in_batch * (GrB_Index) incoming_edges ;
    //     GrB_Index outgoing_edges_in_batch = (GrB_Index) nodes_in_batch * (GrB_Index) outgoing_edges ;
    //     GrB_Index edges_in_batch =  incoming_edges_in_batch + outgoing_edges_in_batch ;
    //     if(edges_in_batch <= 0) {
    //         break ;
    //     }

    //     //---------Outgoing edges section---------
    //     GrB_Index batch_begin = curr_edge ;
    //     GrB_Index batch_end = curr_edge + outgoing_edges_in_batch - 1 ;

    //     GrB_Index range_out[3] ;
    //     range_out[GxB_BEGIN] = batch_begin ;
    //     range_out[GxB_END] = batch_end ;
    //     range_out[GxB_INC] = 1 ;

    //     GRB_TRY (GrB_Vector_clear (Batch_random_out)) ;
    //     GRB_TRY (GrB_Vector_extract (Batch_random_out, NULL, NULL, State, range_out, GxB_RANGE, NULL)) ;
        
    //     Scale_Params scale_out_params ;
    //     scale_out_params.base_edge = (uint64_t) curr_edge ;
    //     GRB_TRY (GrB_Scalar_setElement_UDT (scale_scalar, &scale_out_params)) ;

    //     GRB_TRY (GrB_Vector_clear (Scaled_state)) ;
    //     GRB_TRY (GrB_apply (Scaled_state, NULL, NULL, scale_op, Batch_random_out, scale_scalar, NULL)) ;

    //     GRB_TRY (GrB_Vector_clear (Batch_gather)) ;
    //     GRB_TRY (GrB_Vector_extract (Batch_gather, NULL, NULL, Output_j, Scaled_state, outgoing_edges_in_batch, NULL)) ;

    //     GRB_TRY (GxB_Vector_assign (Output_j, NULL, NULL, Batch_gather, range_out, GxB_RANGE, NULL)) ;
        
    //     Chunk_Index_Params chunk_out_params ;
    //     chunk_out_params.base_node = curr_node ;
    //     chunk_out_params.edges = (uint64_t) outgoing_edges ;
    //     GRB_TRY (GrB_Scalar_setElement_UDT (chunk_index_scalar, &chunk_out_params)) ;

    //     GRB_TRY (GrB_Vector_clear (Batch_chunk)) ;
    //     GRB_TRY (GrB_apply (Batch_chunk, NULL, NULL, chunk_index_op, Batch_random_out, chunk_index_scalar, NULL)) ;

    //     GRB_TRY (GxB_Vector_assign (Output_i, NULL, NULL, Batch_chunk, range_out, GxB_RANGE, NULL)) ;

    //     //---------Incoming edges section---------
    //     if(incoming_edges > 0){
    //         GrB_Index batch_begin = curr_edge + outgoing_edges_in_batch ;
    //         GrB_Index batch_end = curr_edge + edges_in_batch - 1 ;

    //         GrB_Index range_in[3] ;
    //         range_in[GxB_BEGIN] = batch_begin ;
    //         range_in[GxB_END] = batch_end ;
    //         range_in[GxB_INC] = 1 ;

    //         GRB_TRY (GrB_Vector_clear (Batch_random_in)) ;
    //         GRB_TRY (GrB_Vector_extract (Batch_random_in, NULL, NULL, State, range_in, GxB_RANGE, NULL)) ;
            
    //         Scale_Params scale_in_params ;
    //         scale_in_params.base_edge = (uint64_t) curr_edge ;
    //         GRB_TRY (GrB_Scalar_setElement_UDT (scale_scalar, &scale_in_params)) ;

    //         GRB_TRY (GrB_Vector_clear (Scaled_state)) ;
    //         GRB_TRY (GrB_apply (Scaled_state, NULL, NULL, scale_op, Batch_random_in, scale_scalar, NULL)) ;

    //         GRB_TRY (GrB_Vector_clear (Batch_gather)) ;
    //         GRB_TRY (GrB_Vector_extract (Batch_gather, NULL, NULL, Output_i, Scaled_state, incoming_edges_in_batch, NULL)) ;

    //         GRB_TRY (GxB_Vector_assign (Output_i, NULL, NULL, Batch_gather, range_in, GxB_RANGE, NULL)) ;
            
    //         Chunk_Index_Params chunk_out_params ;
    //         chunk_out_params.base_node = curr_node ;
    //         chunk_out_params.edges = (uint64_t) incoming_edges ;
    //         GRB_TRY (GrB_Scalar_setElement_UDT (chunk_index_scalar, &chunk_out_params)) ;

    //         GRB_TRY (GrB_Vector_clear (Batch_chunk)) ;
    //         GRB_TRY (GrB_apply (Batch_chunk, NULL, NULL, chunk_index_op, Batch_random_in, chunk_index_scalar, NULL)) ;

    //         GRB_TRY (GxB_Vector_assign (Output_j, NULL, NULL, Batch_chunk, range_in, GxB_RANGE, NULL)) ;
    //     }

    //     //TODO: section for corners

    //     curr_node += nodes_in_batch ;
    //     curr_edge += edges_in_batch ;
    //     curr_batch_size *= batch_growth_factor ;
    // }

    // //TODO: mirror for undirected case

    // //--------------------------------------------------------------------------
    // // Create scalar for unweighted edges
    // //--------------------------------------------------------------------------

    // GRB_TRY (GrB_Scalar_new(&Scalar_one, GrB_UINT8)) ;
    // GRB_TRY (GrB_Scalar_setElement_UINT8(Scalar_one, 1)) ;

    // //--------------------------------------------------------------------------
    // // Build output matrix
    // //--------------------------------------------------------------------------

    // GRB_TRY (GrB_Matrix_new(&Y, GrB_UINT8, num_nodes, num_nodes)) ;
    // GRB_TRY (GxB_Matrix_build_Scalar_Vector(Y, Output_i, Output_j, Scalar_one, NULL)) ;

    LG_FREE_WORK ;
    (*Yhandle) = Y ;
    return (GrB_SUCCESS) ;
}
