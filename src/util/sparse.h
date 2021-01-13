/**
 * @file   sparse.h
 * @brief  Sparse matrix format definition
 */

#ifndef SPARSE_H_
#define SPARSE_H_


#include "common.h"


namespace stokesdt {

namespace detail {
    
/** 
 * @struct  SparseMatrix
 * @brief   Stores sparse matrix stored in blocked CSR format.
 */
typedef struct SparseMatrix {
    /// the dimension of the block 
    size_t sizeb;
    /// <code>sizeb^2</code>
    size_t sizeb2;
    /// the number of nonzero blocks
    size_t nnzb;
    /// the allocated storage capacity
    size_t maxnnzb;
    /// the number of blocks in row dimension
    size_t nrowb;
    /// the array of nonzero values
    double *val;
    /// the array of block column indices
    size_t *colbidx;
    /// the array of block row pointer
    size_t *rowbptr;
} SparseMatrix;


/// Creates a sparse matrix
bool CreateSparseMatrix(const size_t nrowb,
                        const size_t maxnnzb,
                        const size_t sizeb,
                        SparseMatrix **p_spmat);

/// Resizes the sparse matrix to contain <code>maxnnzb</code> nonzero blocks
bool ResizeSparseMatrix(const size_t maxnnzb,
                        SparseMatrix *spmat);

/// Destroys the sparse matrix
void DestroySparseMatrix(SparseMatrix *spmat);

/// Computes sparse matrix-vector multiplication
void SpMV3x3(const SparseMatrix *spmat,
             const int nrhs,
             const double alpha,
             const int ldx,
             const double *x,
             const double beta,
             const int ldy,
             double *y);

/// Transforms a spase matrix into dense format
void SparseToDense(const SparseMatrix *spmat, const int ldm, double *mat);

} // namespace detail

} // namespace stokesdt


#endif // SPARSE_H_
