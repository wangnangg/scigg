import numpy
import scipy.io
import scipy.linalg
import scipy.sparse.linalg
import time



K = scipy.io.mmio.mmread('spmatrices/circuit_1/circuit_1.mtx').tocsr()
f = numpy.squeeze(scipy.io.mmio.mmread('spmatrices/circuit_1/circuit_1_b.mtx'))

# Direct solve
t0 = time.time()
u = scipy.sparse.linalg.spsolve(K, f)
t1 = time.time()
err = scipy.linalg.norm(K*u - f)
print("spsolve() time = %f seconds" % (t1-t0))
print("err = %e" % err)

t0 = time.time()
u, info = scipy.sparse.linalg.gmres(K, f, restart=30, maxiter=1000, tol=1e-10)
t1 = time.time()
err = scipy.linalg.norm(K*u - f)
print("gmres() time = %f seconds" % (t1-t0))
print("err = %e" % err)
print("info = %d" % info)
