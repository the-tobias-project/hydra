

# cython: profile=False
# copile with py3 compiler.py build_ext --inplace
#changed from https://github.com/pandas-dev/pandas/blob/c1068d9d242c22cb2199156f6fb82eb5759178ae/pandas/_libs/algos.pyx

cimport cython
from cython cimport Py_ssize_t

import numpy as np
cimport numpy as np
from numpy cimport (ndarray,
                    NPY_INT64, NPY_UINT64, NPY_INT32, NPY_INT16, NPY_INT8,
                    NPY_FLOAT32, NPY_FLOAT64,
                    NPY_OBJECT,
                    int8_t, int16_t, int32_t, int64_t, uint8_t, uint16_t,
                    uint32_t, uint64_t, float32_t, float64_t,
                    double_t)

cdef float64_t FP_ERR = 1e-13

cdef extern from "../../../L-BFGS-B-C/src/lbfgsb.h":
#    ctypedef logical
    long setulb(long *n, long *m, double *x,
          double *l, double *u, long *nbd, double *f, double 
          *g, double *factr, double *pgtol, double *wa, long *
          iwa, long *task, long *iprint, long *csave, long *lsave, 
          long *isave, double *dsave)

@cython.boundscheck(False)
@cython.wraparound(False)
def minimize_lbfgs(ndarray[float64_t, ndim=2, mode="c"] C not None, ndarray[float64_t, ndim=1] u not None, 
    ndarray[float64_t, ndim=1]z not None, float64_t rho, ndarray[float64_t, ndim=1, mode="c"] x not None, uint64_t n):

    cdef:
        long max_iter, m, intN, output
        double[:] dsave, low_bnd, upper_bnd, g, wa
        long[:] iwa, isave, csave
        double f, pgtol, factr
        long task
    intN = (<long> n) 
    max_iter = 5000
    #nbd = np.zeros(n, dtype=np.int64)
    cdef long[::1] nbd = np.zeros(n, dtype=np.int64)
    low_bnd = np.zeros(n, dtype=np.float64)
    upper_bnd = np.zeros(n, dtype=np.float64)
    f = 0.0
    #g = np.zeros((n,), np.float64)
    g = np.zeros(n, np.float64)
    wa = np.zeros(25*n + 1180, np.float64) # m is set to 10
    iwa = np.zeros(3*n ,np.int64)
    csave = np.zeros(60,np.int64)
    #cdef cnp.ndarray array = np.array([False,False,False,False], dtype=bool)
    #cdef cnp.uint8_t[:] lsave = np.frombuffer(array, dtype=np.uint8)
    #cdef cnp.ndarray[cnp.uint8_t, ndim = 1, cast=True] lsave
    #lsave=np.array([False, False, False, False], dtype=bool)
    #lsave = np.zeros(4, np.bool)
    cdef long[:] lsave = np.zeros(4, np.int64)
    isave = np.zeros(44, np.int64)
    dsave = np.zeros(29, np.float64)
    task = 1 # START

    m = 10
    output = -1
    pgtol = 1e-5
    factr = 1e7
    n_iterations = 0
    cdef Py_ssize_t i
    
    for n_iterations in range(max_iter):
        setulb(&intN, &m, &x[0],&low_bnd[0], &upper_bnd[0], &nbd[0], &f, &g[0], &factr, &pgtol, &wa[0], &iwa[0], &task,
          &output, &csave[0], &lsave[0], &isave[0], &dsave[0])
        if task >= 10 and task <= 15: #FG range
            v = x + u - z
            eCx = np.exp(np.dot(C,x))
            rhov = rho * v
            f = np.sum(np.log1p(eCx)) + 0.5 * np.sum(rhov*v)
            g = np.dot(C.T, eCx / (1+ eCx)) + rhov
#        elif task == 2:
#            n_iterations += 1
        elif task != 2: 
            return x



    

