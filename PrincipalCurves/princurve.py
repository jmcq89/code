import numpy as np

def kern(x):
    "Gaussian Kernel Profile"
    kx = np.exp(-x/2.0)
    return kx

def p(data, points, h):
    "Evaluate KDE on a set of points based on data"
    data = np.atleast_2d(data)
    n, N = data.shape
    points = np.atleast_2d(points)
    m, M = points.shape
    if m == 1 and M == n: # row vector
        points = np.reshape(points, (n, 1))
        m, M = points.shape
    const = (1.0/N)*((h)**(-n))*(2.0*np.pi)**(-n/2.0)
    probs = np.zeros((M,),dtype = np.float)
    for i in range(M):
        diff = (data - points[:, i, np.newaxis])/h
        x = np.sum(diff*diff, axis = 0)
        probs[i] = np.sum(kern(x), axis = 0)*const
    return probs

def u(data, x, h):
    data = np.atleast_2d(data)
    n, N = data.shape
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = np.reshape(x, (n, 1))
    us = (x - data)/(h**2)
    return us

def c(data, x, h):
    data = np.atleast_2d(data)
    n, N = data.shape
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = x.T
    us = u(data, x, h)
    const = (h**(-n))*(2.0*np.pi)**(-n/2.0)
    u2 = np.sum(us*(us*(h**2)), axis=0)
    cs = kern(u2)*const
    return cs

def p2(data, points, h):
    data = np.atleast_2d(data)
    n, N = data.shape
    points = np.atleast_2d(points)
    m, M = points.shape
    if m == 1 and M == n: # row vector
        points = np.reshape(points, (n, 1))
        m, M = points.shape
    probs = np.zeros((M,),dtype = np.float)
    for i in range(M):
        cs = c(data, points[:, i], h)
        probs[i] = np.sum(cs)
    probs = probs/N
    return probs

def ms(data, x, h):
    "Calculate the mean-shift at point x"
    data = np.atleast_2d(data)
    n, N = data.shape
    const = (1.0/N)*((h)**(-n))*(2.0*np.pi)**(-n/2.0)
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = np.reshape(x, (n, 1))
    unprobs = p(data, x, h)/const
    diff = (data - x)/h
    diff2 = np.sum(diff*diff, axis = 0)
    mx = np.sum(kern(diff2)*data, axis = 1)/unprobs
    mx = np.reshape(mx, (n, 1))
    return mx

def ms2(data, x, h):
    data = np.atleast_2d(data)
    n, N = data.shape
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = np.reshape(x, (n, 1))
    cs = c(data, x, h)
    denom = np.zeros((n, n), dtype = np.float)
    num = np.zeros((n,1), dtype = np.float)
    for i in range(N):
        Sigmainv = np.identity(n)*(1/(h**2))
        denom = denom + cs[i]*Sigmainv
        num = num + cs[i]*(np.dot(Sigmainv, np.reshape(data[:, i], (n, 1))))
    ms = np.dot(np.linalg.inv(denom), num)
    ms = np.reshape(ms, (n, 1))
    return ms

def grad(data, x, h):
    "calculate the local gradient of the kernel density"
    "g(x) = (p(x)/h^2)*[mean-shift(x) - x]"
    data = np.atleast_2d(data)
    n, N = data.shape
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = np.reshape(x, (n, 1))
    probs = p(data, x, h)
    mx = np.reshape(ms2(data, x, h), (n, 1))
    gx = (probs/(h**2))*(mx-x)
    return gx

def grad2(data, x, h):
    "calculate the local gradient of the kernel density"
    "g(x) = (-1/N)sum(ci*ui)"
    data = np.atleast_2d(data)
    n, N = data.shape
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = np.reshape(x, (n, 1))
    cs = c(data, x, h)
    us = u(data, x, h)
    gx = -np.sum(cs*us, axis = 1)/N
    gx = np.reshape(gx, (n, 1))
    return gx

def var(data, x, h):
    "calulates the local ''variance'' at point x"
    "similar to mean-shift in spirit"
    data = np.atleast_2d(data)
    n, N = data.shape
    const = (1.0/N)*((1.0/h)**(n))*(1.0/(2.0*np.pi))**(n/2.0)
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = np.reshape(x, (n, 1))
    unprobs = p(data, x, h)/const
    vs = np.zeros((n, n), dtype = np.float)
    diff = (data - x)/h
    wts = kern(np.sum(diff*diff, axis = 0))
    for i in range(N):
        op = np.dot(np.reshape(diff[:, i]*h, (n, 1)), \
            np.reshape(diff[:, i]*h, (1, n)))
        vs = vs + op*wts[i]
    vs = vs/unprobs
    return vs

def hess(data, x, h):
    "Calculate the local hessian of the kernel density at x"
    "H(x) = (p(x)/h^4)*(var(x) - h^2I_n)"
    data = np.atleast_2d(data)
    n, N = data.shape
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = x.T
    probs = p(data, x, h)
    vs = var(data, x, h)
    In = np.identity(n)
    Hx = (probs/(h**4))*(vs - (h**2)*In)
    return Hx

def hess2(data, x, h):
    data = np.atleast_2d(data)
    n, N = data.shape
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = x.T
    Sigmainv = np.identity(n)*(1/(h**2))
    cs = c(data, x, h)
    us = u(data, x, h)
    Hx = np.zeros((n, n), dtype = np.float)
    for i in range(N):
        op = np.dot(np.reshape(us[:,i], (n, 1)), np.reshape(us[:,i], (1, n)))
        hx = cs[i]*(op - Sigmainv)
        Hx = Hx + hx
    Hx = Hx/N
    return Hx

def localcov(data, x, h):
    "calculate the local (inverse) covariance matrix"
    "Cov(X) = - (1/p(x))*H(X) + (1/p(x)^2)*g(x)*g(x)^T"
    data = np.atleast_2d(data)
    n, N = data.shape
    x = np.atleast_2d(x)
    m, M = x.shape
    if M != 1:
        x = x.T
    probs = p(data, x, h)
    Hx = hess(data, x, h)
    gx = grad(data, x, h)
    Cx = (1/(probs**2))*(np.dot(gx, gx.T)) - (1/probs)*Hx
    return Cx

def princurve(data, points, tol, d, h):
    data = np.atleast_2d(data)
    n, N = data.shape
    points = np.atleast_2d(points)
    m, M = points.shape
    if m == 1 and M == n: # row vector
        points = np.reshape(points, (n, 1))
        m, M = points.shape
    for k in range(M):
        #print "%d of %d" % (k+1, M)
        # if k == int(N/4):
            # print "25% complete"
        # if k == int(N/2):
            # print "50% complete"
        # if k == int(3*N/4):
            # print "75% complete"
        points[:, k] = projectpoint(data, points[:, k], tol, d, h)
    return points
    
def projectpoint(data, point, tol, d, h):
    n, N = data.shape
    converged = False
    iter = 0
    while not(converged):
        H = hess(data, point,h)
        w, v = np.linalg.eigh(H)
        index = np.argsort(w) # arguments that sort from small to large
        V = np.reshape(v[:, index[range(n-d)]], (n, d))
        ospace = np.dot(V, V.T)
        proj = np.reshape(ms(data, point, h), (n, 1)) - np.reshape(point, (n, 1))
        proj = np.dot(ospace, proj) + np.reshape(point, (n, 1))
        diff = np.linalg.norm(np.reshape(point, (n, 1)) - proj)
        point = np.reshape(proj, (n, ))
        iter = iter +1
        if diff < tol:
            converged = True
        if iter > 500:
            converged = True
            print "maximum iterations exceeded"
    return point 