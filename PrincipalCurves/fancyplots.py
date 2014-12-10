from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import princurve as pc
from scipy.spatial.distance import pdist, squareform

script, dset = argv

def distmat(d):
    n, m = d.shape
    if n == 2:
        d = d.T
        n, m = d.shape
    dmat = pdist(d, 'euclidean')
    dmat = squareform(dmat)
    return dmat

def findh(data):
    n, N = data.shape
    hs = np.linspace(N**(-2), 1, num = 20)
    tlprobs = np.zeros(len(hs))
    # for each h calculate the total Leave-one-out log-prob
    for i in range(len(hs)):
        h = hs[i]
        probs = np.zeros(N)
        for t in range(N):
            # probability of data point given data excluding itself
            probs[t] = pc.p(data[:, np.delete(range(N), t)], data[:, t], h)
        tlprobs[i] = np.sum(np.log(probs))
    max = np.argmax(tlprobs)
    h = hs[max]
    return h


def findorder(dmat):
    n, n = dmat.shape
    vals = np.zeros((n, 1), dtype = int)
    out = []
    for i in range(n-1):
        qp = vals[i] # start from last point
        out.extend(qp) # list of variables not allowed to be considered
        drow = dmat[qp, :] 
        drow[:,out] = np.inf # set already seen values, and self, to infinite
        mindist = np.argmin(drow) # smallest remaining distance
        qp = mindist # next point
        vals[i+1] = qp
    return vals

data = np.loadtxt(dset); data = data.T # data set

xmin = data[0, :].min(); xmax = data[0, :].max()
ymin = data[1, :].min(); ymax = data[1, :].max()

points = np.copy(data) # so that manipulations don't change data
points = np.concatenate([points, points+0.01, points - 0.01, points+0.02, points - 0.02, points+0.03, points - 0.03], axis = 1)
tol = (10.0)**(-4); d = 1
n, N = data.shape
h = findh(data)

curve = pc.princurve(data, points, tol, d, h) # O&E Principal Curve
curveorder = findorder(distmat(curve))


f, ax = plt.subplots()
ax.plot(data[0, :], data[1, :], 'kx')
ax.plot(curve[0, curveorder], curve[1, curveorder], 'b.')
plt.show()
