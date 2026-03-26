//------------------------------------------------------------------------------
// LAGraph_RMAT: an implementation of the R-MAT graph generator
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

#define LG_FREE_WORK                        \
{                                           \
    /* free any workspace used here */      \
    GrB_free (&State) ;                     \
    GrB_IndexUnaryOp_free (&RMAT_op_i) ;    \
    GrB_IndexUnaryOp_free (&RMAT_op_j) ;    \
    GrB_Type_free (&param_type) ;           \
    GrB_Scalar_free (&params_scalar) ;      \
    GrB_Type_free (&Uint_Size) ;            \
    GrB_Vector_free (&Output_i) ;           \
    GrB_Vector_free (&Output_j) ;           \
    GrB_Scalar_free (&Scalar_one) ;         \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free (&Y) ;                         \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"

typedef struct
{
    uint64_t a;       // RMAT sector a weight
    uint64_t b;       // RMAT sector b weight
    uint64_t c;       // RMAT sector c weight
    uint64_t d;       // RMAT sector d weight
    uint64_t n;       // Number of bits in address (log_2(num_nodes))
} RMAT_Params;
#define RMAT_PARAMS                                                             \
"typedef struct {                                                            \n"\
"   uint64_t a; uint64_t b; uint64_t c; uint64_t d; uint64_t n;              \n"\
"} RMAT_Params;"

void RMAT_Operator_i_32
(
    void *z,
    const void *x,
    GrB_Index i,
    GrB_Index j,
    const void *y
)
{
    const uint64_t *r = (const uint64_t *) x;
    const RMAT_Params *p = (const RMAT_Params *) y;

    uint64_t state = *r ;
    uint32_t result = 0 ;
    uint64_t ab_threshold = p->a + p->b ;
    
    for(int it = 0; it < p->n; ++it){
        bool cond = (state > ab_threshold) ;
        result = (result << 1) | (cond ? 1 : 0) ;

        // XORshift64 for random state on next iteration
        state ^= state << 13 ;
        state ^= state >> 7 ;
        state ^= state << 17 ;
    }
    *((uint32_t *)z) = result;
}
#define RMAT_I_32                                                         \
"void RMAT_Operator_i_32(                                              \n"\
"    void *z,                                                          \n"\
"    const void *x,                                                    \n"\
"    GrB_Index i,                                                      \n"\
"    GrB_Index j,                                                      \n"\
"    const void *y                                                     \n"\
")                                                                     \n"\
"{                                                                     \n"\
"    const uint64_t *r = (const uint64_t *) x;                         \n"\
"    const RMAT_Params *p = (const RMAT_Params *) y;                   \n"\
"                                                                      \n"\
"    uint64_t state = *r;                                              \n"\
"    uint32_t result = 0;                                              \n"\
"    uint64_t ab_threshold = p->a + p->b;                              \n"\
"                                                                      \n"\
"    for(int it = 0; it < p->n; ++it){                                    \n"\
"        bool cond = (state > ab_threshold);                           \n"\
"        result = (result << 1) | (cond ? 1 : 0);                      \n"\
"        state ^= state << 13;                                         \n"\
"        state ^= state >> 7;                                          \n"\
"        state ^= state << 17;                                         \n"\
"    }                                                                 \n"\
"    *((uint32_t *)z) = result;                                        \n"\
"}                                                                     "

void RMAT_Operator_i_64
(
    void *z,
    const void *x,
    GrB_Index i,
    GrB_Index j,
    const void *y
)
{
    const uint64_t *r = (const uint64_t *) x;
    const RMAT_Params *p = (const RMAT_Params *) y;

    uint64_t state = *r ;
    uint64_t result = 0 ;
    uint64_t ab_threshold = p->a + p->b ;
    
    for(int it = 0; it < p->n; ++it){
        bool cond = (state > ab_threshold) ;
        result = (result << 1) | (cond ? 1 : 0) ;

        // XORshift64 for random state on next iteration
        state ^= state << 13 ;
        state ^= state >> 7 ;
        state ^= state << 17 ;
    }
    *((uint64_t *)z) = result;
}
#define RMAT_I_64                                                         \
"void RMAT_Operator_i_64(                                              \n"\
"    void *z,                                                          \n"\
"    const void *x,                                                    \n"\
"    GrB_Index i,                                                      \n"\
"    GrB_Index j,                                                      \n"\
"    const void *y                                                     \n"\
")                                                                     \n"\
"{                                                                     \n"\
"    const uint64_t *r = (const uint64_t *) x;                         \n"\
"    const RMAT_Params *p = (const RMAT_Params *) y;                   \n"\
"                                                                      \n"\
"    uint64_t state = *r;                                              \n"\
"    uint64_t result = 0;                                              \n"\
"    uint64_t ab_threshold = p->a + p->b;                              \n"\
"                                                                      \n"\
"    for(int it = 0; it < p->n; ++it){                                    \n"\
"        bool cond = (state > ab_threshold);                           \n"\
"        result = (result << 1) | (cond ? 1 : 0);                      \n"\
"        state ^= state << 13;                                         \n"\
"        state ^= state >> 7;                                          \n"\
"        state ^= state << 17;                                         \n"\
"    }                                                                 \n"\
"    *((uint64_t *)z) = result;                                        \n"\
"}                                                                     "

void RMAT_Operator_j_32
(
    void *z,
    const void *x,
    GrB_Index i,
    GrB_Index j,
    const void *y
)
{
    const uint64_t *r = (const uint64_t *) x;
    const RMAT_Params *p = (const RMAT_Params *) y;

    uint64_t state = *r ;
    uint32_t result = 0 ;
    uint64_t ab_threshold = p->a + p->b ;
    uint64_t abc_threshold = ab_threshold + p->c ;
    
    for(int it = 0; it < p->n; ++it){
        bool cond = (
            (state > p->a && state <= ab_threshold) ||
            (state > abc_threshold)
        ) ;
        result = (result << 1) | (cond ? 1 : 0) ;

        // XORshift64 for random state on next iteration
        state ^= state << 13 ;
        state ^= state >> 7 ;
        state ^= state << 17 ;
    }
    *((uint32_t *)z) = result;
}
#define RMAT_J_32                                                           \
"void RMAT_Operator_j_32(                                                \n"\
"    void *z,                                                            \n"\
"    const void *x,                                                      \n"\
"    GrB_Index i,                                                        \n"\
"    GrB_Index j,                                                        \n"\
"    const void *y                                                       \n"\
")                                                                       \n"\
"{                                                                       \n"\
"    const uint64_t *r = (const uint64_t *) x;                           \n"\
"    const RMAT_Params *p = (const RMAT_Params *) y;                     \n"\
"                                                                        \n"\
"    uint64_t state = *r;                                                \n"\
"    uint32_t result = 0;                                                \n"\
"    uint64_t ab_threshold = p->a + p->b;                                \n"\
"    uint64_t abc_threshold = ab_threshold + p->c;                       \n"\
"                                                                        \n"\
"    for(int it = 0; it < p->n; ++it){                                      \n"\
"        bool cond = (                                                   \n"\
"            (state > p->a && state <= ab_threshold) ||                  \n"\
"            (state > abc_threshold)                                     \n"\
"        );                                                              \n"\
"        result = (result << 1) | (cond ? 1 : 0);                        \n"\
"        state ^= state << 13;                                           \n"\
"        state ^= state >> 7;                                            \n"\
"        state ^= state << 17;                                           \n"\
"    }                                                                   \n"\
"    *((uint32_t *)z) = result;                                          \n"\
"}                                                                       "

void RMAT_Operator_j_64
(
    void *z,
    const void *x,
    GrB_Index i,
    GrB_Index j,
    const void *y
)
{
    const uint64_t *r = (const uint64_t *) x;
    const RMAT_Params *p = (const RMAT_Params *) y;

    uint64_t state = *r ;
    uint64_t result = 0 ;
    uint64_t ab_threshold = p->a + p->b ;
    uint64_t abc_threshold = ab_threshold + p->c ;
    
    for(int it = 0; it < p->n; ++it){
        bool cond = (
            (state > p->a && state <= ab_threshold) ||
            (state > abc_threshold)
        ) ;
        result = (result << 1) | (cond ? 1 : 0) ;

        // XORshift64 for random state on next iteration
        state ^= state << 13 ;
        state ^= state >> 7 ;
        state ^= state << 17 ;
    }
    *((uint64_t *)z) = result;
}
#define RMAT_J_64                                                           \
"void RMAT_Operator_j_64(                                                \n"\
"    void *z,                                                            \n"\
"    const void *x,                                                      \n"\
"    GrB_Index i,                                                        \n"\
"    GrB_Index j,                                                        \n"\
"    const void *y                                                       \n"\
")                                                                       \n"\
"{                                                                       \n"\
"    const uint64_t *r = (const uint64_t *) x;                           \n"\
"    const RMAT_Params *p = (const RMAT_Params *) y;                     \n"\
"                                                                        \n"\
"    uint64_t state = *r;                                                \n"\
"    uint64_t result = 0;                                                \n"\
"    uint64_t ab_threshold = p->a + p->b;                                \n"\
"    uint64_t abc_threshold = ab_threshold + p->c;                       \n"\
"                                                                        \n"\
"    for(int it = 0; it < p->n; ++it){                                      \n"\
"        bool cond = (                                                   \n"\
"            (state > p->a && state <= ab_threshold) ||                  \n"\
"            (state > abc_threshold)                                     \n"\
"        );                                                              \n"\
"        result = (result << 1) | (cond ? 1 : 0);                        \n"\
"        state ^= state << 13;                                           \n"\
"        state ^= state >> 7;                                            \n"\
"        state ^= state << 17;                                           \n"\
"    }                                                                   \n"\
"    *((uint64_t *)z) = result;                                          \n"\
"}                                                                       "

int LAGraph_RMAT
(
    // output
    GrB_Matrix *Yhandle,    // Y, created on output
    // input:
    int log2_nodes,
    GrB_Index num_edges, //approximate number of edges
    uint64_t seed,
    double a, double b, double c, double d, // should add to 1.0
    char *msg
)
{

    //--------------------------------------------------------------------------
    // Set up workspace
    //--------------------------------------------------------------------------

    GrB_Vector State = NULL ;
    GrB_IndexUnaryOp RMAT_op_i = NULL ;
    GrB_IndexUnaryOp RMAT_op_j = NULL ;
    GrB_Type param_type = NULL ;
    GrB_Scalar params_scalar = NULL ;
    GrB_Type Uint_Size = NULL ;
    GrB_Vector Output_i = NULL ;
    GrB_Vector Output_j = NULL ;
    GrB_Scalar Scalar_one = NULL ;
    GrB_Matrix Y = NULL ;
    
    LG_CLEAR_MSG ;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_ASSERT (Yhandle != NULL, GrB_NULL_POINTER) ;
    (*Yhandle) = NULL ;
 
    LG_ASSERT_MSG (log2_nodes > 0, GrB_INVALID_VALUE, "log2_nodes must be positive") ;
    LG_ASSERT_MSG (log2_nodes <= 60, GrB_INVALID_VALUE, "log2_nodes cannot exceed 60") ;
    LG_ASSERT_MSG (num_edges > 0, GrB_INVALID_VALUE, "num_edges must be positive") ;
    double sum = a + b + c + d ;
    double tol = 1.0 / (double) ((uint64_t) 1 << 8);
    //FIXME: double check precision and error tolerance
    LG_ASSERT_MSG (fabs(sum - 1.0) <= tol, GrB_INVALID_VALUE, "a,b,c,d must add up to 1.0 (within tolerance 1e-8)") ;
    //--------------------------------------------------------------------------
    // Create random state vector
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&State, GrB_UINT64, num_edges)) ;
    GRB_TRY (GrB_assign (State, NULL, NULL, 0, GrB_ALL, num_edges, NULL)) ;
    LG_TRY (LAGraph_Random_Seed (State, seed, msg)) ;

    //--------------------------------------------------------------------------
    // Create RMAT_Params struct
    //--------------------------------------------------------------------------

    RMAT_Params params = {0} ;
    params.a = (uint64_t) (a * (double) UINT64_MAX) ;
    params.b = (uint64_t) (b * (double) UINT64_MAX) ;
    params.c = (uint64_t) (c * (double) UINT64_MAX) ;
    params.d = (uint64_t) (d * (double) UINT64_MAX) ;
    params.n = log2_nodes ;

    GRB_TRY (GrB_Type_new (&param_type, sizeof(RMAT_Params))) ;
    GRB_TRY (GrB_Scalar_new (&params_scalar, param_type)) ;
    GRB_TRY (GrB_Scalar_setElement_UDT (params_scalar, &params)) ;

    //--------------------------------------------------------------------------
    // Create operators
    //--------------------------------------------------------------------------

    bool use_64 = (log2_nodes > 31) ;
    Uint_Size = (use_64 ? GrB_UINT64 : GrB_UINT32) ;

    if(use_64){
        GRB_TRY (GrB_IndexUnaryOp_new (&RMAT_op_i, (GxB_index_unary_function) RMAT_Operator_i_64, Uint_Size, GrB_UINT64, param_type)) ;
        GRB_TRY (GrB_IndexUnaryOp_new (&RMAT_op_j, (GxB_index_unary_function) RMAT_Operator_j_64, Uint_Size, GrB_UINT64, param_type)) ;
    } else {
        GRB_TRY (GrB_IndexUnaryOp_new (&RMAT_op_i, (GxB_index_unary_function) RMAT_Operator_i_32, Uint_Size, GrB_UINT64, param_type)) ;
        GRB_TRY (GrB_IndexUnaryOp_new (&RMAT_op_j, (GxB_index_unary_function) RMAT_Operator_j_32, Uint_Size, GrB_UINT64, param_type)) ;
    }

    //--------------------------------------------------------------------------
    // Apply operators to I and J vectors
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Vector_new (&Output_i, Uint_Size, num_edges)) ;
    GRB_TRY (GrB_apply (Output_i, NULL, NULL, RMAT_op_i, State, params_scalar, NULL)) ;

    GRB_TRY (GrB_Vector_new (&Output_j, Uint_Size, num_edges)) ;
    GRB_TRY (GrB_apply (Output_j, NULL, NULL, RMAT_op_j, State, params_scalar, NULL)) ;

    //--------------------------------------------------------------------------
    // Create scalar
    //--------------------------------------------------------------------------

    GRB_TRY (GrB_Scalar_new(&Scalar_one, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8(Scalar_one, 1)) ;

    //--------------------------------------------------------------------------
    // Build output matrix
    //--------------------------------------------------------------------------

    GrB_Index num_nodes = ((GrB_Index) 1) << log2_nodes ;
    GRB_TRY (GrB_Matrix_new(&Y, GrB_UINT8, num_nodes, num_nodes)) ;
    GRB_TRY (GxB_Matrix_build_Scalar_Vector(Y, Output_i, Output_j, Scalar_one, NULL)) ;

    LG_FREE_WORK ;
    (*Yhandle) = Y ;
    return (GrB_SUCCESS) ;
}
