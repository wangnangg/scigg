import numpy as np
import scipy.io
import scipy.linalg
import scipy.sparse.linalg
import scipy.optimize
import time
import sys


def test_gmres():
    K = scipy.io.mmio.mmread(sys.argv[1]).tocsc()
    f = np.squeeze(scipy.io.mmio.mmread(sys.argv[2]))
    M2 = scipy.sparse.linalg.spilu(K)
    M_x = lambda x: M2.solve(x)
    M = scipy.sparse.linalg.LinearOperator(K.shape, M_x)
    t0 = time.time()
    u, info = scipy.sparse.linalg.gmres(K, f, restart=int(sys.argv[3]), M=M, maxiter=int(sys.argv[4]), tol=1e-10)
    t1 = time.time()
    err = scipy.linalg.norm(K*u - f)
    print("gmres() time = %f seconds" % (t1-t0))
    print("err = %e" % err)
    print("info = %d" % info)

def test_min_bgfs():
    def obj_fun(x):
        print("trying: ", x)
        return x[0] * x[0] + 2.0 * (x[1] - 1) * (x[1] - 1) + 5.0
    def obj_der(x):
        grad = np.zeros_like(x)
        grad[0] = 2 * x[0]
        grad[1] = 4 * (x[1] - 1)
        return grad
    x0 = np.array([10, 10])
    res = scipy.optimize.minimize(obj_fun, x0, method='BFGS', jac=obj_der, 
            options={'disp': True, 'gtol': 1e-10})
    print(res)

test_min_bgfs()
