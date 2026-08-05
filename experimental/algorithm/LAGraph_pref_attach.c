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
//TODO: document

#define LG_FREE_WORK                            \
{                                               \
    GrB_free (&State) ;                         \
    GrB_free (&State2) ;                        \
    GrB_Vector_free (&Scaled_state) ;           \
    GrB_Vector_free (&Output_i) ;               \
    GrB_Vector_free (&Output_j) ;               \
    GrB_Scalar_free (&Scalar_one) ;             \
    GrB_Vector_free (&Batch_random) ;           \
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
    GrB_free (&M) ;                         \
    LAGraph_Delete (G, msg) ;               \
}

#include "LG_internal.h"
#include "LAGraphX.h"

typedef struct
{
    uint64_t base_edge ;
    uint64_t base_node ;
    uint64_t edges_per_node ;
} Scale_Params ;
#define SCALE_PARAMS                \
"typedef struct {                   \n"\
"    uint64_t base_edge ;           \n"\
"    uint64_t base_node ;           \n"\
"    uint64_t edges_per_node ;      \n"\
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

// TODO: comment
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
    uint64_t node_offset = (uint64_t) i / p->edges_per_node ;
    uint64_t edge_count = p->base_edge + node_offset * p->edges_per_node ;
    if (edge_count == 0){
        *((uint64_t *) z) = 0 ;
    } else {
        *((uint64_t *) z) = (*r) % edge_count ;
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
"    uint64_t node_offset = (uint64_t) i / p->edges_per_node ;          \n"\
"    uint64_t edge_count = p->base_edge + node_offset * p->edges_per_node ; \n"\
"    if (edge_count == 0){                                              \n"\
"        *((uint64_t *) z) = 0 ;                                        \n"\
"    } else {                                                            \n"\
"        double scaled = ((double)(*r) / (double) UINT64_MAX)    \n"\
"            * (double)(edge_count - 1) ;                                \n"\
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

//LAGr_pref_attach
//  takes seed matrix in form of GrB_Matrix
//LAGraph_pref_attach
//  creates seed
//
int LAGraph_pref_attach
(
    // output
    LAGraph_Graph *G,    // G, created on output
    // input:
    GrB_Index num_nodes,
    uint64_t seed,
    GrB_Index batch_size,
    double batch_growth_factor,
    GrB_Vector Input_i,
    GrB_Vector Input_j,
    GrB_Index incoming_edges,
    GrB_Index outgoing_edges,
    // FIXME: is this valirable name and comment incorrect?
    // TODO: change to enum
    bool directed, // if directed, do outgoing edges and mirror
    char *msg
)
{

    //--------------------------------------------------------------------------
    // Set up workspace
    //--------------------------------------------------------------------------

    GrB_Vector State = NULL ;
    GrB_Vector State2 = NULL ;
    GrB_Vector Scaled_state = NULL ;
    GrB_Vector Output_i = NULL ;
    GrB_Vector Output_j = NULL ;
    GrB_Scalar Scalar_one = NULL ;
    GrB_Vector Batch_random = NULL ;
    GrB_Vector Batch_gather = NULL ;
    GrB_Vector Batch_chunk = NULL ;
    GrB_IndexUnaryOp scale_op = NULL ;
    GrB_IndexUnaryOp chunk_index_op = NULL ;
    GrB_Type scale_param_type = NULL ;
    GrB_Type chunk_index_param_type = NULL ;
    GrB_Scalar scale_scalar = NULL ;
    GrB_Scalar chunk_index_scalar = NULL ;
    GrB_Matrix M = NULL ;

    LG_CLEAR_MSG ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_ASSERT (G != NULL, GrB_NULL_POINTER) ;
    (*G) = NULL ;

    LG_ASSERT (Input_i != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (Input_j != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT_MSG (num_nodes > 0, GrB_INVALID_VALUE, "num_nodes must be positive") ;
    LG_ASSERT_MSG (outgoing_edges > 0, GrB_INVALID_VALUE, "outgoing_edges must be positive") ;
    LG_ASSERT_MSG (incoming_edges >= 0, GrB_INVALID_VALUE, "incoming_edges must be non-negative") ;
    LG_ASSERT_MSG (batch_size > 0, GrB_INVALID_VALUE, "batch_size must be positive") ;
    LG_ASSERT_MSG (batch_growth_factor >= 1, GrB_INVALID_VALUE, "batch_growth_factor must be >= 1") ;

    //--------------------------------------------------------------------------
    // Evaluate seed vectors
    //--------------------------------------------------------------------------

    // FIXME: seed is not a good name for this
    GrB_Index seed_nvals_i ;
    GRB_TRY (GrB_Vector_nvals (&seed_nvals_i, Input_i)) ;
    GrB_Index seed_nvals_j ;
    GRB_TRY (GrB_Vector_nvals (&seed_nvals_j, Input_j)) ;
    LG_ASSERT_MSG (seed_nvals_i == seed_nvals_j, GrB_DIMENSION_MISMATCH, "Input vectors must have same number of edges") ;

    uint64_t max_i = 0 ;
    uint64_t max_j = 0 ;
    GRB_TRY (GrB_reduce (&max_i, NULL, GrB_MAX_MONOID_UINT64, Input_i, NULL)) ;
    GRB_TRY (GrB_reduce (&max_j, NULL, GrB_MAX_MONOID_UINT64, Input_j, NULL)) ;
    uint64_t base_node = ((max_i > max_j) ? max_i : max_j) + 1 ;

    //--------------------------------------------------------------------------
    // Calculate num_edges and num_batches
    //--------------------------------------------------------------------------

    if(!directed){
        incoming_edges = 0 ;
    }

    uint64_t num_nodes_to_add = (num_nodes < base_node ? 0 : num_nodes - base_node) ;
    uint64_t num_edges = seed_nvals_i ;

    uint64_t num_batches = 0 ;
    uint64_t nodes_counted = 0 ;
    uint64_t b_size = batch_size ;

    while(nodes_counted < num_nodes_to_add){
        uint64_t nodes_in_batch = b_size ;

        if (nodes_counted + nodes_in_batch > num_nodes_to_add){
            nodes_in_batch = num_nodes_to_add - nodes_counted ;
        }

        num_edges += nodes_in_batch * (incoming_edges + outgoing_edges) ;
        nodes_counted += nodes_in_batch ;
        ++num_batches ;
        b_size = (b_size * batch_growth_factor == 0 ? 1 : b_size * batch_growth_factor) ;
        // printf ("%" PRIu64 " counted out of %" PRIu64 "\n", (uint64_t) nodes_counted, (uint64_t) num_nodes_to_add) ;
    }

    //--------------------------------------------------------------------------
    // Initialize Output_i and Output_j vectors
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&Output_i, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&Output_j, GrB_UINT64, num_edges)) ;

    GrB_Index range_seed[3] ;
    range_seed[GxB_BEGIN] = 0 ;
    range_seed[GxB_END] = seed_nvals_i - 1 ;
    range_seed[GxB_INC] = 1 ;

    GRB_TRY (GrB_Vector_assign (Output_i, NULL, NULL, Input_i, range_seed, GxB_RANGE, NULL)) ;
    GRB_TRY (GrB_Vector_assign (Output_j, NULL, NULL, Input_j, range_seed, GxB_RANGE, NULL)) ;

    //--------------------------------------------------------------------------
    // Create random state vector
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&State, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_assign (State, NULL, NULL, 0, GrB_ALL, num_edges, NULL)) ;
    LG_TRY (LAGraph_Random_Seed (State, seed, msg)) ;

    //--------------------------------------------------------------------------
    // Prepare parameter structs, operators, and workspace vectors
    //--------------------------------------------------------------------------

    // GrB_set (GrB_GLOBAL, true, GxB_BURBLE) ;
    GRB_TRY (GxB_Type_new (
        &scale_param_type, sizeof (Scale_Params), "Scale_Params", SCALE_PARAMS
    )) ;
    // TODO: use GxB method to feed in the JIT strings
    // verify JIT is running with the burble
    // looking for JIT compile then JIT cpu NOT generic
    GRB_TRY (GrB_Type_new (&chunk_index_param_type, sizeof (Chunk_Index_Params))) ;

    GRB_TRY (GrB_Scalar_new (&scale_scalar, scale_param_type)) ;
    GRB_TRY (GrB_Scalar_new (&chunk_index_scalar, chunk_index_param_type)) ;

    GRB_TRY (GrB_IndexUnaryOp_new (&scale_op, (GxB_index_unary_function) pref_attach_scale_op, GrB_UINT64, GrB_UINT64, scale_param_type)) ;
    GRB_TRY (GrB_IndexUnaryOp_new (&chunk_index_op, (GxB_index_unary_function) pref_attach_chunk_index_op, GrB_UINT64, GrB_UINT64, chunk_index_param_type)) ;

    // TODO: maybe rename vectors to lowercase
    GRB_TRY (GrB_Vector_new (&Scaled_state, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&Batch_random, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&Batch_gather, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_Vector_new (&Batch_chunk, GrB_UINT64, num_edges)) ;

    //--------------------------------------------------------------------------
    // Begin main loop
    //--------------------------------------------------------------------------

    uint64_t curr_node = base_node ;
    GrB_Index curr_edge = seed_nvals_i ;
    uint64_t curr_batch_size = batch_size ;

    for(int batch = 0; batch < num_batches; ++batch){
        // Prepare variables for batch
        uint64_t nodes_in_batch = LAGRAPH_MIN (num_nodes - curr_node, curr_batch_size) ;
        GrB_Index outgoing_edges_in_batch = (GrB_Index) nodes_in_batch * (GrB_Index) outgoing_edges ;
        GrB_Index incoming_edges_in_batch = (GrB_Index) nodes_in_batch * (GrB_Index) incoming_edges ;
        GrB_Index edges_in_batch = incoming_edges_in_batch + outgoing_edges_in_batch ;
        if(edges_in_batch <= 0) {
            break ;
        }

        //---------Outgoing edges section---------
        if(outgoing_edges_in_batch > 0){
            // TODO: consider removing Batch_random, Scaled_state, Batch_chunk
            // May require some funky arithmetic but it should be all linear
            GRB_TRY (GrB_Vector_resize (Batch_random, outgoing_edges_in_batch)) ;
            GRB_TRY (GrB_Vector_resize (Scaled_state, outgoing_edges_in_batch)) ;
            GRB_TRY (GrB_Vector_resize (Batch_gather, outgoing_edges_in_batch)) ;
            GRB_TRY (GrB_Vector_resize (Batch_chunk, outgoing_edges_in_batch)) ;

            GrB_Index batch_begin = curr_edge ;
            GrB_Index batch_end = curr_edge + outgoing_edges_in_batch - 1 ;

            GrB_Index range_out[3] ;
            range_out[GxB_BEGIN] = batch_begin ;
            range_out[GxB_END] = batch_end ;
            range_out[GxB_INC] = 1 ;

            Chunk_Index_Params chunk_out_params ;
            chunk_out_params.base_node = curr_node ;
            chunk_out_params.edges = (uint64_t) outgoing_edges ;
            GRB_TRY (GrB_Scalar_setElement_UDT (chunk_index_scalar, &chunk_out_params)) ;

            GRB_TRY (GrB_assign (Batch_chunk, NULL, NULL, (uint64_t) 0, GrB_ALL, outgoing_edges_in_batch, NULL)) ;
            GRB_TRY (GrB_apply (Batch_chunk, NULL, NULL, chunk_index_op, Batch_chunk, chunk_index_scalar, NULL)) ;

            GRB_TRY (GrB_Vector_assign (Output_i, NULL, NULL, Batch_chunk, range_out, GxB_RANGE, NULL)) ;
            GRB_TRY (GrB_Vector_assign (Output_j, NULL, NULL, Batch_chunk, range_out, GxB_RANGE, NULL)) ;

            GRB_TRY (GrB_Vector_clear (Batch_random)) ;
            GRB_TRY (GrB_Vector_extract (Batch_random, NULL, NULL, State, range_out, GxB_RANGE, NULL)) ;

            Scale_Params scale_out_params ;
            scale_out_params.base_edge = (uint64_t) curr_edge ;
            scale_out_params.base_node = curr_node ;
            scale_out_params.edges_per_node = (uint64_t) outgoing_edges ;
            GRB_TRY (GrB_Scalar_setElement_UDT (scale_scalar, &scale_out_params)) ;

            // scale down random numbers in about range(0, i)
            GRB_TRY (GrB_apply (Batch_random, NULL, NULL, scale_op,
                Batch_random, scale_scalar, NULL)) ;

            GxB_print(Batch_gather, GxB_SUMMARY) ;
            GxB_print(Output_j, GxB_SUMMARY) ;
            GxB_print(Scaled_state, GxB_SUMMARY) ;
            GRB_TRY (GxB_Vector_extract_Vector (Batch_gather, NULL, NULL, Output_j, Scaled_state, NULL)) ;
            GRB_TRY (GrB_Vector_assign (Output_j, NULL, NULL, Batch_gather, range_out, GxB_RANGE, NULL)) ;

            GRB_TRY (GrB_Vector_clear (Batch_gather)) ;
        }

        //---------Incoming edges section---------
        if(incoming_edges_in_batch > 0){
            GRB_TRY (GrB_Vector_resize (Batch_random, incoming_edges_in_batch)) ;
            GRB_TRY (GrB_Vector_resize (Scaled_state, incoming_edges_in_batch)) ;
            GRB_TRY (GrB_Vector_resize (Batch_gather, incoming_edges_in_batch)) ;
            GRB_TRY (GrB_Vector_resize (Batch_chunk, incoming_edges_in_batch)) ;

            GrB_Index batch_begin = curr_edge + outgoing_edges_in_batch ;
            GrB_Index batch_end = batch_begin + incoming_edges_in_batch - 1 ;

            GrB_Index range_in[3] ;
            range_in[GxB_BEGIN] = batch_begin ;
            range_in[GxB_END] = batch_end ;
            range_in[GxB_INC] = 1 ;

            Chunk_Index_Params chunk_in_params ;
            chunk_in_params.base_node = curr_node ;
            chunk_in_params.edges = (uint64_t) incoming_edges ;
            GRB_TRY (GrB_Scalar_setElement_UDT (chunk_index_scalar, &chunk_in_params)) ;

            GRB_TRY (GrB_assign (Batch_chunk, NULL, NULL, (uint64_t) 0, GrB_ALL, incoming_edges_in_batch, NULL)) ;
            GRB_TRY (GrB_apply (Batch_chunk, NULL, NULL, chunk_index_op, Batch_chunk, chunk_index_scalar, NULL)) ;

            GRB_TRY (GrB_Vector_assign (Output_i, NULL, NULL, Batch_chunk, range_in, GxB_RANGE, NULL)) ;
            GRB_TRY (GrB_Vector_assign (Output_j, NULL, NULL, Batch_chunk, range_in, GxB_RANGE, NULL)) ;

            GRB_TRY (GrB_Vector_clear (Batch_random)) ;
            GRB_TRY (GrB_Vector_extract (Batch_random, NULL, NULL, State, range_in, GxB_RANGE, NULL)) ;

            Scale_Params scale_in_params ;
            scale_in_params.base_edge = (uint64_t) (curr_edge + outgoing_edges_in_batch) ;
            scale_in_params.base_node = curr_node ;
            scale_in_params.edges_per_node = (uint64_t) incoming_edges ;
            GRB_TRY (GrB_Scalar_setElement_UDT (scale_scalar, &scale_in_params)) ;

            GRB_TRY (GrB_Vector_clear (Scaled_state)) ;
            GRB_TRY (GrB_apply (Scaled_state, NULL, NULL, scale_op, Batch_random, scale_scalar, NULL)) ;

            GRB_TRY (GxB_Vector_extract_Vector (Batch_gather, NULL, NULL, Output_i, Scaled_state, NULL)) ;
            GRB_TRY (GrB_Vector_assign (Output_i, NULL, NULL, Batch_gather, range_in, GxB_RANGE, NULL)) ;
            GRB_TRY (GrB_Vector_clear (Batch_gather)) ;
        }

        curr_node += nodes_in_batch ;
        curr_edge += edges_in_batch ;
        curr_batch_size *= batch_growth_factor ;
    }

    //--------------------------------------------------------------------------
    // Mirror edges if undirected
    //-------------------------------------------------------------------------

    if(!directed){
        GrB_Index total_edges ;
        GRB_TRY (GrB_Vector_nvals (&total_edges, Output_i));

        GRB_TRY (GrB_Vector_resize (Batch_gather, total_edges)) ;
        GRB_TRY (GrB_Vector_resize (Batch_chunk, total_edges)) ;

        GrB_Index range_source[3] ;
        range_source[GxB_BEGIN] = 0 ;
        range_source[GxB_END] = total_edges - 1 ;
        range_source[GxB_INC] = 1 ;

        GRB_TRY (GrB_Vector_clear (Batch_gather)) ;
        GRB_TRY (GrB_Vector_extract (Batch_gather, NULL, NULL, Output_j, range_source, GxB_RANGE, NULL)) ;

        GRB_TRY (GrB_Vector_clear (Batch_chunk)) ;
        GRB_TRY (GrB_Vector_extract (Batch_chunk, NULL, NULL, Output_i, range_source, GxB_RANGE, NULL)) ;

        GRB_TRY (GrB_Vector_resize (Output_i, total_edges * 2)) ;
        GRB_TRY (GrB_Vector_resize (Output_j, total_edges * 2)) ;

        GrB_Index range_dest[3] ;
        range_dest[GxB_BEGIN] = total_edges ;
        range_dest[GxB_END] = (total_edges * 2) - 1 ;
        range_dest[GxB_INC] = 1 ;

        GRB_TRY (GrB_Vector_assign (Output_i, NULL, NULL, Batch_gather, range_dest, GxB_RANGE, NULL)) ;
        GRB_TRY (GrB_Vector_assign (Output_j, NULL, NULL, Batch_chunk, range_dest, GxB_RANGE, NULL)) ;
    }

    //--------------------------------------------------------------------------
    // Create scalar for unweighted edges
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Scalar_new(&Scalar_one, GrB_BOOL)) ;
    GRB_TRY (GrB_Scalar_setElement_BOOL (Scalar_one, true)) ;

    //--------------------------------------------------------------------------
    // Build output Graph
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Matrix_new(&M, GrB_BOOL, num_nodes, num_nodes)) ;
    GRB_TRY (GxB_Matrix_build_Scalar_Vector(M, Output_i, Output_j, Scalar_one, NULL)) ;

    LAGraph_Kind kind = directed ? LAGraph_ADJACENCY_DIRECTED : LAGraph_ADJACENCY_UNDIRECTED ;
    LG_TRY (LAGraph_New (G, &M, kind, msg)) ;
    M = NULL ;
    // TODO: possibly redundant?
    LG_TRY (LAGraph_DeleteSelfEdges (*G, msg)) ;

    LG_FREE_WORK ;

    return (GrB_SUCCESS) ;
}
