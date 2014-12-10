from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import princurve as pc
from scipy.spatial.distance import pdist, squareform
import os

script, dset, sax = argv

def distmat(d):
    n, m = d.shape
    if n == 2:
        d = d.T
        n, m = d.shape
    dmat = pdist(d, 'euclidean')
    dmat = squareform(dmat)
    return dmat

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

script_dir = os.path.dirname(__file__)
rel_path = "CurveData/" + dset 
abs_file_path = os.path.join(script_dir, rel_path)
os.chdir(abs_file_path) # CD to location of Data 

# Load in Data Set
data = np.loadtxt('sample.dta'); data = data.T # data set
g = np.loadtxt('gencurve.dta'); g = g.T # generating curve
hs = np.loadtxt('hscurve.dta'); hs = hs.T # Hastie-Stuetzle Principal Curve
pg = np.loadtxt('project.dta'); pg = pg.T # Polygonal Alg. Principal Curve

# CD back to princurve location 
os.chdir('../..')

xmin = data[0, :].min(); xmax = data[0, :].max()
ymin = data[1, :].min(); ymax = data[1, :].max()

pg = pg.T
pg = pg[pg[:,sax].argsort()]
pg = pg.T
pgorder = findorder(distmat(pg))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(data[0, :], data[1, :], 'kx')
ax.plot(hs[0, :], hs[1, :], 'y')
ax.plot(pg[0, :], pg[1, :], 'b.')
plt.axis([xmin, xmax, ymin, ymax])
plt.show()

## Find the O&E Principal Curve
points = np.copy(data) # so that manipulations don't change data
tol = (10.0)**(-4); d = 1
n, N = data.shape
h1 = N**(-1./(n)) # Scott's Rule: n+4
h2 = N**(-1./(n+1))
h3 = N**(-1./(n+1.75))

curve1 = pc.princurve(data, points, tol, d, h1) # O&E Principal Curve
points = np.copy(data) # so that manipulations don't change data
curve2 = pc.princurve(data, points, tol, d, h2) # O&E Principal Curve
points = np.copy(data) # so that manipulations don't change data
curve3 = pc.princurve(data, points, tol, d, h3) # O&E Principal Curve

curve1 = curve1.T; curve1 = curve1[curve1[:, sax].argsort()]; curve1 = curve1.T
curve1order = findorder(distmat(curve1))
curve2 = curve2.T; curve2 = curve2[curve2[:, sax].argsort()]; curve2 = curve2.T
curve2order = findorder(distmat(curve2))
curve3 = curve3.T; curve3 = curve3[curve3[:, sax].argsort()]; curve3 = curve3.T
curve3order = findorder(distmat(curve3))


X, Y = np.mgrid[xmin:xmax:500j, ymin:ymax:500j]
points = np.vstack([X.ravel(), Y.ravel()])
n, N = data.shape

kernel = pc.p(data, points, h1)
Z1 = np.reshape(kernel.T, X.shape)
kernel = pc.p(data, points, h2)
Z2 = np.reshape(kernel.T, X.shape)
kernel = pc.p(data, points, h3)
Z3 = np.reshape(kernel.T, X.shape)

#, extent=[xmin, xmax, ymin, ymax]
f, (ax1,ax2, ax3) = plt.subplots(1, 3)
ax1.imshow(np.rot90(Z1), cmap=plt.cm.Greens, extent=[xmin, xmax, ymin, ymax])
ax1.plot(data[0, :], data[1, :], 'kx')
ax1.plot(curve1[0, curve1order], curve1[1, curve1order], 'b.')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, ymax)

ax2.imshow(np.rot90(Z2), cmap=plt.cm.Greens, extent=[xmin, xmax, ymin, ymax])
ax2.plot(data[0, :], data[1, :], 'kx')
ax2.plot(curve2[0, curve2order], curve2[1, curve2order], 'b.')
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(ymin, ymax)

ax3.imshow(np.rot90(Z3), cmap=plt.cm.Greens, extent=[xmin, xmax, ymin, ymax])
ax3.plot(data[0, :], data[1, :], 'kx')
ax3.plot(curve3[0, curve3order], curve3[1, curve3order], 'b.')
ax3.set_xlim(xmin, xmax)
ax3.set_ylim(ymin, ymax)
plt.show()
