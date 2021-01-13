/**
 * @file   sparse.cc
 * @brief  Sparse matrix operations
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <x86intrin.h>
#include "sparse.h"
#include "log.h"


namespace stokesdt {

namespace detail {

bool CreateSparseMatrix(const size_t nrowb,
                        const size_t maxnnzb,
                        const size_t sizeb,
                        SparseMatrix **p_spmat)
{
    SparseMatrix *spmat = new SparseMatrix;
    spmat->maxnnzb = maxnnzb;
    spmat->nnzb = 0;
    spmat->nrowb = nrowb;
    spmat->sizeb = sizeb;
    spmat->sizeb2 = sizeb * sizeb;
    spmat->rowbptr = (size_t *)AlignMalloc(sizeof(size_t) * nrowb);
    spmat->colbidx = (size_t *)malloc(sizeof(size_t) * maxnnzb);
    spmat->val = (double *)malloc(sizeof(double) * maxnnzb * sizeb * sizeb);
    if (NULL == spmat->rowbptr ||
        NULL == spmat->colbidx ||
        NULL == spmat->val) {
        LOG_ERROR("Failed to allocate memory: %lld.\n",
            sizeof(double) * maxnnzb * sizeb * sizeb);
        return false;
    }

    *p_spmat = spmat;

    return true;
}


bool ResizeSparseMatrix(const size_t maxnnzb,
                        SparseMatrix *spmat)
{
    if (spmat->maxnnzb < maxnnzb) {
        size_t sizeb = spmat->sizeb;
        spmat->maxnnzb = maxnnzb;           
        spmat->colbidx = (size_t *)realloc(spmat->colbidx, sizeof(size_t) * maxnnzb);
        spmat->val = (double *)malloc(sizeof(double) * maxnnzb * sizeb * sizeb);
        if (NULL == spmat->colbidx ||
            NULL == spmat->val) {
            LOG_ERROR("Failed to allocate memory: %lld.\n",
                sizeof(double) * maxnnzb * sizeb * sizeb);
            return false;
        }
    }
    
    return true;
}


void DestroySparseMatrix(SparseMatrix *spmat)
{
    if (spmat != NULL) {
        AlignFree(spmat->rowbptr);
        free(spmat->colbidx);
        free(spmat->val);
        delete spmat;
    }
}


void SpMV3x3(const SparseMatrix *spmat,
             const int nrhs,
             const double alpha,
             const int ldx,
             const double *x,
             const double beta,
             const int ldy,
             double *y)
{
    int nrowb = spmat->nrowb;
    size_t *rowbptr = spmat->rowbptr;
    size_t *colbidx = spmat->colbidx;
    double *val = spmat->val;   
    #pragma omp parallel for
    for (size_t i = 0; i < nrowb; i++) {
        __declspec(align(kAlignLen))
            double tmp[kAlignLen/sizeof(double)];        
        size_t startb = rowbptr[i];
        size_t endb = rowbptr[i + 1];
        size_t row = 3 * i;
        double _beta = beta;
        /* compute a row */
        for (size_t j = startb; j < endb; j++) {
            size_t col = 3 * colbidx[j];       
        #if defined(__SSE3__)
            // SSE3 kernel
            __m128d a0010 =_mm_loadu_pd (&(val[j * 9 + 0]));
            __m128d a0111 =_mm_loadu_pd (&(val[j * 9 + 3]));
            __m128d a0212 =_mm_loadu_pd (&(val[j * 9 + 6]));
            __m128d a20 =_mm_set1_pd (val[j * 9 + 2]);
            __m128d a21 =_mm_set1_pd (val[j * 9 + 5]);
            __m128d a22 =_mm_set1_pd (val[j * 9 + 8]);
            for (size_t k = 0; k < nrhs; k++) {
                __m128d vx0 =_mm_set1_pd (x[k * ldx + col + 0]);
                __m128d vx1 =_mm_set1_pd (x[k * ldx + col + 1]);
                __m128d vx2 =_mm_set1_pd (x[k * ldx + col + 2]);
                __m128d vy01 = _mm_mul_pd (a0010, vx0);
                __m128d vy2  = _mm_mul_pd (a20, vx0);
                vy01 = _mm_add_pd (vy01, _mm_mul_pd (a0111, vx1));
                vy2  = _mm_add_pd (vy2,  _mm_mul_pd (a21, vx1));
                vy01 = _mm_add_pd (vy01, _mm_mul_pd (a0212, vx2));
                vy2  = _mm_add_pd (vy2,  _mm_mul_pd (a22, vx2));
                _mm_store_pd(&(tmp[0]), vy01);
                _mm_store_pd(&(tmp[2]), vy2);
                y[k * ldy + row + 0] = 
                    _beta * y[k * ldy + row + 0] + alpha * tmp[0];
                y[k * ldy + row + 1] =
                    _beta * y[k * ldy + row + 1] + alpha * tmp[1];
                y[k * ldy + row + 2] =
                    _beta * y[k * ldy + row + 2] + alpha * tmp[2];
            }
        #else
            // scalar kernel    
            for (size_t k = 0; k < nrhs; k++) {
                tmp[0]  = val[j * 9 + 0] * x[k * ldx + col + 0];
                tmp[1]  = val[j * 9 + 1] * x[k * ldx + col + 0];
                tmp[2]  = val[j * 9 + 2] * x[k * ldx + col + 0];
                tmp[0] += val[j * 9 + 3] * x[k * ldx + col + 1];
                tmp[1] += val[j * 9 + 4] * x[k * ldx + col + 1];
                tmp[2] += val[j * 9 + 5] * x[k * ldx + col + 1];
                tmp[0] += val[j * 9 + 6] * x[k * ldx + col + 2];
                tmp[1] += val[j * 9 + 7] * x[k * ldx + col + 2];
                tmp[2] += val[j * 9 + 8] * x[k * ldx + col + 2];
                y[k * ldy + row + 0]
                    = _beta * y[k * ldy + row + 0] + alpha * tmp[0];
                y[k * ldy + row + 1] 
                    = _beta * y[k * ldy + row + 1] + alpha * tmp[1];
                y[k * ldy + row + 2] 
                    = _beta * y[k * ldy + row + 2] + alpha * tmp[2];
            }
        #endif          
            _beta = 1.0;
        } /* for (j = startb; j < endb; j++) */
    } /* for (i = 0; i < nrowb; i++) */
}


void SparseToDense(const SparseMatrix *spmat, const int ldm, double *mat)
{
    size_t nrowb = spmat->nrowb;
    size_t sizeb = spmat->sizeb;
    size_t sizeb2 = spmat->sizeb2;
    size_t *rowbptr = spmat->rowbptr;
    size_t *colbidx = spmat->colbidx;
    double *val = spmat->val;

    memset (mat, 0, sizeof(double) * ldm * sizeb * nrowb);
    for (size_t i = 0; i < nrowb; i++) {
        size_t start = rowbptr[i];
        size_t end = rowbptr[i + 1];
        for (size_t j = start; j < end; j++) {
            for (size_t p = 0; p < sizeb; p++) {
                size_t row = i * sizeb + p;
                for (size_t q = 0; q < sizeb; q++) {
                    size_t col = colbidx[j] * sizeb + q;
                    size_t idx = p * sizeb + q;
                    mat[row * ldm + col] = val[j * sizeb2 + idx];
                }
            }
        }
    }   
}

} // namespace detail

} // namespace stokesdt
