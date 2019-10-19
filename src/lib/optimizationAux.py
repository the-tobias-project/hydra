# Optimization helper functions

import numpy as np
from sklearn.metrics import log_loss
from numpy import sum,maximum, exp, log, log1p
from scipy.optimize import fmin_l_bfgs_b as bfgs
from scipy.optimize.lbfgsb import _minimize_lbfgsb
from scipy.optimize.optimize import MemoizeJac, wrap_function
from scipy.optimize._lbfgsb import setulb
from numpy.linalg import cholesky
from numpy.core._methods import _sum
from numpy.core import umath
#from itertools import imap,starmap,izip
from operator import mul
from numpy.linalg import _umath_linalg
#from profilehooks import profile
import numba as nb


from scipy.special import expit

import pdb
umr_sum = umath.add.reduce
add     = np.add
sub     = np.subtract
dot     = np.core.multiarray.dot
trans   = np.transpose
div     = np.divide
asarray = np.asarray
zeros   = np.zeros
int32   = np.int32
float64 = np.float64
array   = np.array
mult    = np.multiply
npabs    = np.abs
eye     = np.eye
empty   = np.empty
lstsq   = _umath_linalg.lstsq_m

# ADMM helpers
#@nb.jit(nb.types.Tuple((nb.f8, nb.f8[:,:])) (nb.f8[:], nb.f8[:,:], nb.f8[:,:],nb.f8[:,:],nb.f8, nb.i8), locals={'v':nb.f8[:,:],
#  'eCx': nb.f8[:,:], 'f': nb.f8, 'g': nb.f8[:,:], 'y':nb.f8[:,:]})
def l2_log(x, C, z, u, rho, n):
  """ assumes C is an numpy array"""
  v = sub(add(x,u), z)
  eCx = exp(dot(C,x))

  f = _sum(log1p(eCx)) + 0.5 * rho * _sum(mult(v, v))
  g = dot(C.T, div(eCx, (add(1, eCx)))) +  rho * v
  return f, g

def bfgs_update(C, u, z, rho, x0):
  n = len(x0)
  args = (C, z, u, rho,n)
  # bfgs(func, x0, args, bounds, callback)
  return bfgs(l2_log, x0, args=args, bounds=None, callback=callback)

def bfgs_gutted(C, u, z, rho, x0):
   args = (C, z, u, rho)
   fun = MemoizeJac(l2_log)
   jac = fun.derivative
   return _minimize_lbfgsb(fun, x0, args, jac, None, None, 10,
      2.2204460492503131e-09, 1e-5, 1e-8, 15000, 15000, -1, None, 20)['x']

#@numba.jit(cache=True)
def bfgs_more_gutted(C, u, z, rho, x, n):
  #max_iter = 15000
  max_iter = 3000
  nbd = zeros(n, int32)
  low_bnd = zeros(n, float64)
  upper_bnd = zeros(n, float64)

  x = array(x, float64)
  f = 0.0
  g = zeros((n,), float64)
  wa = zeros(20*n + 5*n + 1180, float64)
  iwa = zeros(3*n, int32)
  task = zeros(1, 'S60')
  csave = zeros(1, 'S60')
  lsave = zeros(4, int32)
  isave = zeros(44, int32)
  dsave = zeros(29, float64)
  task[:] = 'START'

  n_iterations = 0
  while 1:
    setulb(10, x, low_bnd, upper_bnd, nbd, f, g, 1e7,
        1e-5, wa, iwa, task, -1, csave, lsave,isave, dsave, 20)
    task_str = task.tostring()
    if task_str.startswith(b'FG'):
      v = sub(add(x,u),z)
      eCx = exp(dot(C,x))
      rhov = rho * v
      f = umr_sum(log1p(eCx), None, None, None, False) + 0.5 * umr_sum(mul(rhov,v), None, None, None, False)
      g = dot(C.T, div(eCx, (add(1, eCx)))) + rhov
    elif task_str.startswith(b'NE'):
      n_iterations += 1
      if n_iterations >= max_iter:
        return x
    else:
      return x


def simple_newton (C, u, z, rho, x, n):
    max_iter = 20
    alpha = 0.1
    BETA  = 0.5
    TOLERANCE = 1e-5
    I = np.eye(n)
    umz = sub(u,z)

    for i in range(max_iter):
        v = add(x,umz)
        eCx = exp(dot(C,x))
        rhov = rho * v

        f = umr_sum(log1p(eCx), None, None, None, False) + 0.5 * umr_sum(mul(rhov,v), None, None, None, False)
        logitECX = div(eCx, add(1,eCx))
        g = dot(C.T, logitECX) + rhov
        H = dot(C.T,  np.diag(logitECX/(1 + eCx))).dot(C) + rho*I
        dx = np.linalg.lstsq(-H, g)[0]
        dfx = g.T.dot(dx)
        if npabs(dfx) < TOLERANCE:
            break
        t = 1
        while True:
            xnew = x+t*dx
            v = add(xnew,umz)
            eCx = exp(dot(C,xnew))
            fnew = umr_sum(log1p(eCx), None, None, None, False) + 0.5 * rho * umr_sum(mul(v,v), None, None, None, False)
            if fnew > f + alpha*t*dfx:
                t *= BETA
            else:
                break
        x = xnew
    return x

nb.jit(nb.f8[:](nb.f8[:,:], nb.f8[:], nb.f8[:], nb.f8, nb.f8[:], nb.i8),
    locals={'max_iter': nb.i8, 'alpha':nb.f8, 'BETA':nb.f8, 'TOLERANCE': nb.f8, 'umz': nb.f8[:],
      'm': nb.i8, 'f':nb.f8, 't':nb.f8, 'dfx': nb.f8, 'v':nb.f8[:], 'eCx':nb.f8[:],
      'dx': nb.f8[:], 'fnew':nb.f8},
    fastmath={'fast'}, nogil=True,
    cache=True,nopython=True)
def other_newton (C, u, z, rho, x, n):
    max_iter = 5
    alpha = 0.1
    BETA  = 0.5
    TOLERANCE = 1e-5

    umz = sub(u,z)
    m = C.shape[0]
    v = add(x,umz)
    eCx = exp(dot(C,x))
    f = umr_sum(log1p(eCx), None, None, None, False) + rho * 0.5 * umr_sum(mul(v,v), None, None, None, False)
    #f = sum(log1p(eCx)) + rho * 0.5 * sum(mul(v,v))

    for i in range(max_iter):
        logitECX = div(eCx, add(1,eCx))
        g = dot(C.T, logitECX) + rho * v
        g = g[:, np.newaxis]
        #H = dot(dot(C.T,  np.diag(logitECX/(1 + eCx))), C) + rho*I
        H = hess(C,logitECX/(1 + eCx),n,m, rho)
        dx = lstsq(-H, g, 1e-4, signature='ddd->ddid')[0]
        dfx = dot(g.T, dx)
        if npabs(dfx) < TOLERANCE:
            break
        t = 1
        while True:
            xnew = x+t*dx.T
            xnew = xnew[0]
            v = add(xnew,umz)
            eCx = exp(dot(C,xnew))
            fnew = umr_sum(log1p(eCx), None, None, None, False) + 0.5 * rho * umr_sum(mul(v,v), None, None, None, False)
            if fnew > f + alpha*t*dfx:
                t *= BETA
            else:
                break
        x = xnew
        if (f-fnew) < TOLERANCE:
            break
        f = fnew
    return x

@nb.jit(nb.f8[:, :](nb.f8[:,:], nb.f8[:], nb.i8, nb.i8, nb.f8),
    locals={'i': nb.i4, 'l': nb.i4, 'val': nb.f8, 'H': nb.f8[:,:]},
    cache=True, nopython=True)
def hess(C, M, n, m, rho):
    H = empty((n,n))
    for i in range(n):
        for l in range(i):
            val = 0
            #H[i,l] += dot(M, C[:,i]* C[:,l])
            for j in range(m):
                val += M[j]*C[j,i]* C[j,l]
            H[l,i] = val
            H[i,l] = val
    for i in range(n):
        val = 0
        for j in range(m):
            val += M[j]*C[j,i]* C[j,i]
        H[i, i] = val + rho
    return H


@nb.jit(nb.types.Tuple((nb.f8[:,:], nb.f8[:,:], nb.f8[:,:], nb.f8))(nb.f8[:,:], nb.f8[:], nb.i8, nb.i8, nb.f8),
    locals={'i': nb.i4, 'l': nb.i4, 'val': nb.f8, 'H': nb.f8[:,:], 'M': nb.f8[:], 'logitECX': nb.f8[:], 'eCx': nb.f8[:]},
    cache=True, nopython=True, fastmath={'fast'})
def ltri_Hessians(C, x, n, m, rho):
    eCx = exp(dot(C,x))
    f = sum(log1p(eCx)) #+ rho * 0.5 * umr_sum(mul(v,v), None, None, None, False)
    logitECX = div(eCx, add(1,eCx))
    g = dot(C.T, logitECX).reshape((-1,1))# + rho * v
    M=logitECX/(1 + eCx)
    H = zeros((n,n))
    for i in range(n):
        for l in range(i):
            val = 0
            for j in range(m):
                val += M[j]*C[j,i]* C[j,l]
            H[i,l] = val
    d = np.empty((1, n))
    for i in range(n):
        val = 0
        for j in range(m):
            val += M[j]*C[j,i]* C[j,i]
        d[0, i] = val
    return H, d, g.T, f


@nb.jit(nb.f8(nb.f8[:,:], nb.f8[:]))
def function_values(C, x):
    eCx = exp(dot(C,x))
    f = sum(log1p(eCx))
    return f


def shrinkage(a, kappa):
  return maximum(0, a - kappa) - maximum(0, -a - kappa)


def l1_OLS(A, b, lam, x, z):
  return 0.5 * sum((A.dot(x) - b) ** 2 ) + lam * np.norm(z,1)

def lasso_admm_cholesky(A, rho):
  n, p = A.shape
  if n > p:
    chol = cholesky(A.T.dot(A) + rho * np.eye(p, dtype=np.float))
  else:
    chol = cholesky(np.eye(n, dtype=np.float) + 1.0 / rho * A.dot(A.T))
  return chol

def callback(x):
  pass




