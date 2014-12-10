import numpy as np
import matplotlib.pyplot as plt
import princurve as pc
import waveletdenoise as wdn
from sys import argv
import time

script, J, ns, ms = argv

J = int(J); ns = int(ns); ms = int(ms)

# generate the noisy signal
def signal(n, noise):
    "define random points around a sin curve"
    t = np.linspace(0.0, 4*np.pi, num = n)
    D = -np.sin(0.75*t) - 0.5*np.cos(1.5*t) - 0.25*np.sin(5*t)
    e = np.random.normal(scale = noise, size = n)
    x = D + e
    SNR = np.sum(D**2, axis = 0)/(n*noise**2)
    return t, x, SNR, D

# LOOCV max probability to find h
def findh(data):
    n, N = data.shape
    hs = np.linspace(N**(-1), 1, num = 15)
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

def findhthr(data, tol):
    n, N = data.shape
    # start with scott's rule:
    sigma = N**(-1./(n+4))
    print sigma
    denom = N*(N-1)*n
    const = ((sigma)**(-n))*(2.0*np.pi)**(-n/2.0)
    for t in range(10):
        print t
        probs = np.zeros(N)
        for i in range(N):
            # calculate LOOProb:
            prob = pc.p(data[:, np.delete(range(N), i)], data[:, i], sigma)
            # calculate innter term:
            diff = (data[:, np.delete(range(N), i)] - data[:, i])
            dist2 = diff*diff
            Gij = const*np.exp(-dist2/(2*sigma**2))
            inner = np.sum(Gij*dist2)
            probs[i] = inner/prob
        tprobs = np.sum(probs)
        print tprobs
        sigma = tprobs/denom
        print sigma
    return sigma

# Scaling filter for LA(8) wavelet
g = [-0.0757657147893407, -0.0296355276459541,
    0.4976186676324578, 0.8037387518052163,
    0.2978577956055422, -0.0992195435769354,
    -0.0126039672622612, 0.0322231006040713]

N = 2**J
lsigma = np.linspace(-12, 2, num = ns)
sigma = np.exp(lsigma)/(1+np.exp(lsigma))

MSEW = np.zeros((ms, ns))
MSEP = np.zeros((ms, ns))
SNR = np.zeros(ns)

for m in range(ms):
    print 'Monte Carlo Simulation %d of %d' %(m+1, ms)
    for i in range(ns):
        # print 'iteration %d of %d, sigma = %f' %(i+1, ns, sigma[i])
        t, x, SNR[i], D = signal(N, sigma[i])
        hdata = np.vstack([x])
        data = np.vstack([t, x])
        points = np.copy(data) # so that manipulations don't change data

        # wavelet de-nolising
        Xwav = wdn.wden(x, g)

        # principal curve
        # print 'selecting h'
        tol = (10.0)**(-3); d = 1; h = findh(hdata)
        # print 'finding principal curve'
        curve = pc.princurve(data, points, tol, d, h)
        curve = curve.T; curve = curve[curve[:, 0].argsort()]; curve = curve.T

        # MSE from true signal
        MSEW[m, i] = np.sum((Xwav - D)**2)/N
        MSEP[m, i] = np.sum((curve[1, :] - D)**2)/N

lMSEWav = np.mean(np.log(MSEW), axis = 0)
lMSEPav = np.mean(np.log(MSEP), axis = 0)
lMSEWsd = np.std(np.log(MSEW), axis = 0)
lMSEPsd = np.std(np.log(MSEW), axis = 0)
lMSEWub = lMSEWav + 2*lMSEWsd
lMSEWlb = lMSEWav - 2*lMSEWsd
lMSEPub = lMSEPav + 2*lMSEPsd
lMSEPlb = lMSEPav - 2*lMSEPsd

t, x, snr, D = signal(N, 0.1)
Sdb = 10*np.log10(SNR)

f, ax = plt.subplots(2)
ax[0].plot(t, D, 'b.')
ax[0].set_title('True Signal: N = %d' %N)
ax[0].set_xlim(t.min(), t.max())
ax[1].plot(Sdb, lMSEWav, 'b', label = 'Wavelet')
ax[1].plot(Sdb, lMSEWub,'b--')
ax[1].plot(Sdb, lMSEWlb,'b--')
ax[1].plot(Sdb, lMSEPav, 'k', label = 'Principal Curve')
ax[1].plot(Sdb, lMSEPub,'k--')
ax[1].plot(Sdb, lMSEPlb,'k--')
ax[1].set_title('Wavelet vs. Principal Component De-noising')
ax[1].set_ylabel('log-Mean-Squared Error')
ax[1].set_xlabel('SNR (in dB)')
ax[1].set_xlim(Sdb.min(), Sdb.max())
legend = ax[1].legend(loc='upper right')
# plt.show()
pltname = 'denoiseplotN'+str(N)+'.png'
f.savefig(pltname)

