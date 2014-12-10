import numpy as np

def scal2wav(g):
    h = []
    for l in range(len(g)):
        h.extend([(-1)**(l)*g[len(g)-1-l]])
    return h

def pyramidstep(Vj, h, g):
    L = len(h); M = len(Vj);
    W = np.zeros(M/2); V = np.zeros(M/2);
    for t in range(M/2):
        u = 2*t + 1
        W[t] = h[0]*Vj[u]
        V[t] = g[0]*Vj[u]
        for n in range(1, L):
            u = u-1
            if u < 0:
                u = M-1
            W[t] = W[t] + h[n]*Vj[u]
            V[t] = V[t] + g[n]*Vj[u]
    return W, V 

def ipyramidstep(Wj, Vj, h, g):
    M = len(Wj); L = len(h)
    l = -2; m = -1;
    V = np.zeros(2*M)
    for t in range(M):
        l = l+2; m = m+2; u = t; i = 1; k=0;
        V[l] = h[i]*Wj[u]+g[i]*Vj[u]
        V[m] = h[k]*Wj[u]+g[k]*Vj[u]
        if L > 2:
            for n in range(1, L/2):
                u = u +1
                if u >= M:
                    u = 0
                i = i + 2; k = k+2;
                V[l] = V[l] + h[i]*Wj[u]+g[i]*Vj[u]
                V[m] = V[m] + h[k]*Wj[u]+g[k]*Vj[u]
    return V

def dwt(X, h, g):
    N = len(X)
    J = np.log2(N)
    Vj = X; W = np.array([])
    for j in range(int(J)):
        Wj, Vj = pyramidstep(Vj, h, g)
        W = np.concatenate([W, Wj])
    W = np.concatenate([W, Vj])
    return W
    
def idwt(W, h, g):
    N = len(W); J = int(np.log2(N))
    WJ = {}
    e = 0
    for j in range(J):
        s = e
        e = s + N/(2**(j+1))
        WJ[j] = W[range(s, e)]
    Vj = W[range(e,N)]
    for j in range(J):
        Wj = WJ[J - j-1]
        Vj = ipyramidstep(Wj, Vj, h, g)
    X = Vj
    return X

def mad(X, h, g):
    W1, V = pyramidstep(X, h, g)
    sigma = np.median(np.abs(W1))/0.6745
    return sigma

def hthresh(w, sigma, N):
    delta = sigma*np.sqrt(2*np.log(N))
    wt = w if np.abs(w) > delta else 0
    return wt

def wden(X, g):
    h = scal2wav(g)
    N = len(X)
    W = dwt(X, h, g)
    WT = np.zeros(N)
    sigma = mad(X, h, g)
    # don't threshold the scaling coeff:
    for i in range(N-1):
        WT[i] = hthresh(W[i], sigma, N)
    WT[N-1] = W[N-1]
    XH = idwt(WT, h, g)
    return XH