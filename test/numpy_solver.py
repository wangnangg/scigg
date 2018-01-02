import numpy
import scipy.io
import scipy.linalg
import scipy.sparse.linalg
import time
import sys



K = scipy.io.mmio.mmread(sys.argv[1]).tocsc()
f = numpy.squeeze(scipy.io.mmio.mmread(sys.argv[2]))
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
