################################################################
## Python port of https://github.com/markusheinonen/adaptivegp
## (commit ed9ac09, Oct 18 2019) by @dereksoeder
##
## adaptivegp.py
## Last updated July 25, 2024
################################################################

from itertools import permutations
import numpy as np
from numpy import random
from numpy.linalg import cond
from scipy.linalg import cholesky, eigvals, inv, solve, solve_triangular
from scipy.sparse import csr_matrix as sparse
from scipy.spatial.distance import cdist
from scipy.stats import gmean


#//#
#//# MATLAB compatibility functions
#//#

def MLlength(a):
    return max(a.shape)  # note that :mat:`length(a)` is NOT the same as :mat:`size(a, 1)`!

def MLisfield(a, b):
    #return hasattr(a, b)
    return False  # :mat:`isfield` always returns 0 when called on a class (which `gpmodel` is)

def MLslash(a, b):
    if np.isscalar(b):
        return a / b
    return solve(b.T, a.T).T  # `solve(b.T, a.T)` returns `x.T` such that `b.T @ x.T == a.T` (derives from `x == a @ inv(b)` --> `x @ b == a`)

def MLbackslash(a, b):
    if np.isscalar(a):
        return b / a
    return solve(a, b)  # `solve(a, b)` returns `x` such that `a @ x == b` (derives from `x == inv(a) @ b`)

def import_plt(f):
    # function decorator that causes matplotlib.pyplot to be imported on demand, so we can operate without it as long as plotting functionality isn't used
    # (only the functions decorated with `@import_plt` will require pyplot)

    def wrapper(*args, **kwargs):
        global plt
        try:
            assert plt
        except(NameError):
            import matplotlib.pyplot as plt
        return f(*args, **kwargs)

    return wrapper

@import_plt
def MLfigure():
    return plt.figure()

@import_plt
def MLdrawnow():
    plt.draw_all()

@import_plt
def plt_show():
    plt.show()


#//#
#//# computegrads.m
#//#

def computegrads(pars):
    n = MLlength(pars.ell)

    #//# numpy doesn't allow expanding an array by assigning out of bounds, so preallocate instead
    nd = ('l' in pars.nsfuncs)*(n-1) + 1 + \
        ('s' in pars.nsfuncs)*(n-1) + 1 + \
        ('o' in pars.nsfuncs)*(n-1) + 1 + \
        ('m' in pars.nsfuncs)*3 + \
        ('a' in pars.nsfuncs)*3 + \
        ('b' in pars.nsfuncs)*3

#    nd = 1

    # compute derivs
    j = 0
    grad = np.zeros(nd)
    if 'l' in pars.nsfuncs:
        grad[j:j+n] = deriv_ell(pars)
        j += n
    else:
        grad[j] = np.mean(deriv_ell(pars,1))
        j += 1
    if 's' in pars.nsfuncs:
        grad[j:j+n] = deriv_sigma(pars)
        j += n
    else:
        grad[j] = np.mean(deriv_sigma(pars,1))
        j += 1
    if 'o' in pars.nsfuncs:
        grad[j:j+n] = deriv_omega(pars)
        j += n
    else:
        grad[j] = np.mean(deriv_omega(pars,1))
        j += 1
    if 'm' in pars.nsfuncs:
        grad[j:j+3] = deriv_mus(pars)
        j += 3
    if 'a' in pars.nsfuncs:
        grad[j:j+3] = deriv_alphas(pars)
        j += 3
    if 'b' in pars.nsfuncs:
        grad[j:j+3] = deriv_betas(pars)

    return grad


#//#
#//# computemse.m
#//#

def computemse(ypred, yts):
    return np.mean((yts - ypred)**2, axis=0)


#//#
#//# computenlpd.m
#//#

def computenlpd(ypred, yts, var):
    nlpd = np.zeros(MLlength(yts))
    for i in range(MLlength(yts)):
        nlpd[i] = logmvnpdf(yts[np.newaxis,i], ypred[np.newaxis,i], var[np.newaxis,i])[0]
    nlpd = np.mean(nlpd)
#    nlpd = - nlpd / MLlength(yts)

    return nlpd


#//#
#//# computenmse.m
#//#

def computenmse(ypred, yts, ytr):
    return np.mean((yts - ypred)**2, axis=0) / np.mean((yts - np.mean(yts))**2, axis=0)


#//#
#//# denormalise.m
#//#

def denormalise(pars, ft=0, ftstd=0, lt=0, st=0, ot=0): #, ftderiv=0, ftderivstd=0):
    # denormalisation
    xtr = pars.xtr * pars.xscale + pars.xbias
    ytr = pars.ytr * pars.yscale + pars.ybias
    ft = ft * pars.yscale + pars.ybias
    ftstd = ftstd * pars.yscale
    lt = lt * gmean(pars.xscale)
    ot = ot * pars.yscale
    st = st * pars.yscale

#    ftderiv = ftderiv * pars.yscale + pars.ybias
#    ftderivstd = ftderivstd * pars.yscale

    return xtr, ytr, ft, ftstd, lt, st, ot


#//#
#//# deriv_alphas.m
#//#

def deriv_alphas(gp):
# derivative of the alphas over MLL

    Kl = gausskernel(gp.xtr,gp.xtr,gp.betaell,   gp.alphaell, gp.tol)
    Ks = gausskernel(gp.xtr,gp.xtr,gp.betasigma, gp.alphasigma, gp.tol)
    Ko = gausskernel(gp.xtr,gp.xtr,gp.betaomega, gp.alphaomega, gp.tol)

    a_l = MLbackslash(Kl,(gp.l_ell   - gp.l_muell))
    a_s = MLbackslash(Ks,(gp.l_sigma - gp.l_musigma))
    a_o = MLbackslash(Ko,(gp.l_omega - gp.l_muomega))

    dKl = 2 / gp.alphaell   * Kl
    dKs = 2 / gp.alphasigma * Ks
    dKo = 2 / gp.alphaomega * Ko

    dl = 0.5*np.sum(np.diag((a_l @ a_l.T - inv(Kl))*dKl))
    ds = 0.5*np.sum(np.diag((a_s @ a_s.T - inv(Ks))*dKs))
    do = 0.5*np.sum(np.diag((a_o @ a_o.T - inv(Ko))*dKo))

    return dl, ds, do


#//#
#//# deriv_betas.m
#//#

def deriv_betas(gp):
# derivative of the betas over MLL

    Kl = gausskernel(gp.xtr,gp.xtr,gp.betaell,   gp.alphaell, gp.tol)
    Ks = gausskernel(gp.xtr,gp.xtr,gp.betasigma, gp.alphasigma, gp.tol)
    Ko = gausskernel(gp.xtr,gp.xtr,gp.betaomega, gp.alphaomega, gp.tol)

    a_l = MLbackslash(Kl,(gp.l_ell   - gp.l_muell))
    a_s = MLbackslash(Ks,(gp.l_sigma - gp.l_musigma))
    a_o = MLbackslash(Ko,(gp.l_omega - gp.l_muomega))

    dKl = gp.betaell**(-3)   * gp.D * Kl
    dKs = gp.betasigma**(-3) * gp.D * Ks
    dKo = gp.betaomega**(-3) * gp.D * Ko

    dl = 0.5*np.sum(np.diag((a_l @ a_l.T - inv(Kl))*dKl))# + 0.5/pars.betaell - 2 # gamma prior
    ds = 0.5*np.sum(np.diag((a_s @ a_s.T - inv(Ks))*dKs))# + 0.5/pars.betasigma - 2
    do = 0.5*np.sum(np.diag((a_o @ a_o.T - inv(Ko))*dKo))# + 0.5/pars.betaomega - 2

    return dl, ds, do


#//#
#//# deriv_ell.m
#//#

def deriv_ell(pars, scalar=False):
# derivative of the ell latent function wrt MLL

    if ('a' in pars.nsfuncs) or ('b' in pars.nsfuncs):
        pars.Kl = gausskernel(pars.xtr, pars.xtr, pars.betaell, pars.alphaell, pars.tol)

    n = MLlength(pars.xtr)
    Ky = nsgausskernel(pars.xtr,pars.xtr,pars.l_ell, pars.l_ell, pars.l_sigma, pars.l_sigma, pars.l_omega)
    a = MLbackslash(Ky,pars.ytr)
    A = np.outer(a,a) - inv(Ky)

    # correct in the log-transform?
    if MLlength(pars.ell) == 1 or scalar:
        Kf = nsgausskernel(pars.xtr, pars.xtr, pars.l_ell, pars.l_ell, pars.l_sigma, pars.l_sigma, -np.inf)
        dKl = np.exp(np.mean(pars.l_ell))**(-2) * (pars.D * Kf)

        dl_l = 0.5*(np.diag( A@dKl )) - MLbackslash(pars.Kl,(pars.l_ell - pars.l_muell))

        dl_l = np.full(n, np.sum(dl_l))

        if MLisfield(pars, 'Ll'):
            dl_l = MLbackslash(pars.Ll,dl_l)
    else:
        ell = pars.ell
        sigma = pars.sigma

        # compute dK matrix first, then cross-slice from it
        L = ell[:,np.newaxis]**2 + ell[np.newaxis]**2
        E = np.exp(-pars.D/L)
        R = np.sqrt( 2*np.outer(ell,ell) / L )
        dK = np.outer(ell,ell) * np.outer(sigma,sigma) * E / R * (L**(-3)) * (4 * pars.D * ell[:,np.newaxis]**2 - ell[:,np.newaxis]**4 + ell[np.newaxis]**4)

#        dK = np.zeros((n,n))
#
#        for i in range(n):
#            li = ell[i]
#            for j in range(n):
#                lj = ell[j]
#                lij = li**2 + lj**2
#                eij = np.exp(-pars.D[i,j]/lij)
#                rij = np.sqrt( (2*li*lj)/(li**2 + lj**2) )
#                sij = sigma[i]*sigma[j]
#                dK[i,j] = li*lj * sij*eij/rij*lij**(-3) * (4*pars.D[i,j]*li**2 - li**4 + lj**4)

        dl_l = np.zeros(n)
        for i in range(n):
#            ei = np.zeros(n)
#            ei[i] = 1
#            Mi = ei[:,np.newaxis] | ei[np.newaxis]

            # original code:
#            dl[i] = 0.5*np.sum(np.diag(A @ (Mi * dK) ))
            # optimised:
            # i) replace np.diag() with np.sum(,axis=1) since we only need diagonal results
            # ii) make Mi*dK sparse
#           dl[i] = 0.5*np.sum( np.sum(sparse(Mi * dK).T.multiply(A),axis=1) )
            # iii) construct sparse matrix directly from the col/row
            dl_l[i] = 0.5*np.sum( np.sum( sparse((
                                              np.append(dK[:,i], dK[i,:]),
                                              ( np.append(np.arange(n), np.full(n,i)),
                                                np.append(np.full(n,i), np.arange(n)) )
                                          )).T.multiply(A),
                                          axis=1 ) )

        dl_l -= MLbackslash(pars.Kl,(pars.l_ell - pars.l_muell))
#        dl -= MLbackslash(Kl,(np.log(ell)-np.mean(np.log(ell), axis=0)))

        # transform
#        if MLisfield(pars, 'Ll'):
#            dl = pars.Ll.T@dl

    if 'l' in pars.nsfuncs:
        dwl_l = pars.Ls.T@dl_l
    else:
        dwl_l = MLbackslash(pars.Ls,dl_l)

    return dwl_l


#//#
#//# deriv_mus.m
#//#

def deriv_mus(gp):
# derivative of the mean parameters against MLL

    # update if necessary
    if ('a' in gp.nsfuncs) or ('b' in gp.nsfuncs):
        gp.Kl = gausskernel(gp.xtr,gp.xtr,gp.betaell,   gp.alphaell,   gp.tol)
        gp.Ks = gausskernel(gp.xtr,gp.xtr,gp.betasigma, gp.alphasigma, gp.tol)
        gp.Ko = gausskernel(gp.xtr,gp.xtr,gp.betaomega, gp.alphaomega, gp.tol)

    dl = np.sum(MLbackslash(gp.Kl,(gp.l_ell - gp.l_muell)))
    ds = np.sum(MLbackslash(gp.Ks,(gp.l_sigma - gp.l_musigma)))
    do = np.sum(MLbackslash(gp.Ko,(gp.l_omega - gp.l_muomega)))

    return dl, ds, do


#//#
#//# deriv_omega.m
#//#

def deriv_omega(pars, scalar=False):
# derivative of the noise latent function over the MLL
#

    n = MLlength(pars.xtr)

    if ('a' in pars.nsfuncs) or ('b' in pars.nsfuncs):
        pars.Ko = gausskernel(pars.xtr,pars.xtr,pars.betaomega, pars.alphaomega, pars.tol)

    Ky = nsgausskernel(pars.xtr,pars.xtr,pars.l_ell, pars.l_ell, pars.l_sigma, pars.l_sigma, pars.l_omega)
    a = MLbackslash(Ky,pars.ytr)
    A = np.outer(a,a) - inv(Ky)

    dK = np.diag( 2.*np.exp(2*pars.l_omega) )
    dl_o = 0.5*np.diag(A@dK) - MLbackslash(pars.Ko,(pars.l_omega - pars.l_muomega))

    if scalar:
        dl_o = np.full(n, np.sum(dl_o))

#    if MLisfield(pars, 'Lo'):
    if 'o' in pars.nsfuncs:
        dwl_o = pars.Lo.T@dl_o
    else:
        dwl_o = MLbackslash(pars.Lo,dl_o)

    return dwl_o


#//#
#//# deriv_sigma.m
#//#

def deriv_sigma(gp, scalar=False):
# derivative of the sigma latent function wrt MLL

    n = MLlength(gp.xtr)

    if ('a' in gp.nsfuncs) or ('b' in gp.nsfuncs):
        gp.Ks = gausskernel(gp.xtr,gp.xtr,gp.betasigma, gp.alphasigma, gp.tol)

    Ky = nsgausskernel(gp.xtr, gp.xtr, gp.l_ell, gp.l_ell, gp.l_sigma, gp.l_sigma, gp.l_omega)
    Kf = nsgausskernel(gp.xtr, gp.xtr, gp.l_ell, gp.l_ell, gp.l_sigma, gp.l_sigma, -np.inf)

    a = MLbackslash(Ky,gp.ytr)
    A = np.outer(a,a) - inv(Ky)

    dl_s = 2*np.diag( A @ Kf ) - MLbackslash(gp.Ks,(gp.l_sigma - gp.l_musigma))

    if scalar:
        dl_s = np.full(n, np.sum(dl_s))

    if 's' in gp.nsfuncs:
        dwl_s = gp.Ls.T@dl_s
    else:
        dwl_s = MLbackslash(gp.Ls,dl_s)

    return dwl_s


#//#
#//# gausskernel.m
#//#

def gausskernel(X1,X2, ell=1, sigma=1, omega=0):
# non-stationary (scalar) gaussian kernel
    if MLlength(X1) == MLlength(X2):
        K = sigma**2 * np.exp(-0.5* cdist(X1,X2,'sqeuclidean') / ell**2) + omega**2*np.eye(X1.shape[0])
    else:
        K = sigma**2 * np.exp(-0.5* cdist(X1,X2,'sqeuclidean') / ell**2)

    return K


#//#
#//# gpmll.m
#//#

def gpmll(x, y, pars):
    return logmvnpdf(y, np.zeros(MLlength(x)), gausskernel(x,x,np.mean(pars.ell, axis=0), np.mean(pars.sigma, axis=0), np.mean(pars.eps, axis=0)))


#//#
#//# gpmodel.m
#//#

class gpmodel:
    # constructor
    def __init__(self, x, y):
        # kernel parameters
        self.init_ell = 0.05
        self.init_sigma = 0.30
        self.init_omega = 0.05
        self.betaell = 0.20
        self.betasigma = 0.20
        self.betaomega = 0.30
        self.alphaell = 1
        self.alphasigma = 1
        self.alphaomega = 1
        self.muf = 0
        self.tol = 1e-3

        # remaining parameters
        self.nsfuncs = 'lso'
        self.optim = 'grad'
        self.warmups = 300
        self.iters = 250
        self.delta = 0.65
        self.maxdepth = 9
        self.leapfrog = 0.02
        self.restarts = 3
        self.graditers = 5000
        self.plotiters = False
        self.verbose = True

        # dimensions
        self.n = x.shape[0]
        self.d = x.shape[1]
        self.p = y.shape[1]

        # normalise inputs to [0,1] and output to max(abs(y)) = 1
        self.xbias = np.min(x, axis=0)
        self.xscale = np.max(x, axis=0) - np.min(x, axis=0)
        self.yscale = np.max(np.abs(y-np.mean(y, axis=0)), axis=0)
        self.ybias = np.mean(y, axis=0)
        x = (x - self.xbias) / self.xscale
        y = (y - self.ybias) / self.yscale
        self.muf = np.mean(y, axis=0)

        self.xtr = x
        self.ytr = y
        self.D = cdist(self.xtr, self.xtr, 'sqeuclidean')
        self.n = self.xtr.shape[0]

        self.init()

    # init function
    def init(self):
        # prior covariance matrices
        self.Kl = gausskernel(self.xtr, self.xtr, self.betaell,   self.alphaell,   self.tol)
        self.Ks = gausskernel(self.xtr, self.xtr, self.betasigma, self.alphasigma, self.tol)
        self.Ko = gausskernel(self.xtr, self.xtr, self.betaomega, self.alphaomega, self.tol)

        # cholesky decompositions of these
        self.Ll = cholesky(self.Kl, lower=True)
        self.Ls = cholesky(self.Ks, lower=True)
        self.Lo = cholesky(self.Ko, lower=True)

        # set initial parameters in white-log domain
        self.wl_ell = solve_triangular( self.Ll, np.full(self.n, np.log(self.init_ell)), lower=True )
        self.wl_sigma = solve_triangular( self.Ls, np.full(self.n, np.log(self.init_sigma)), lower=True )
        self.wl_omega = solve_triangular( self.Lo, np.full(self.n, np.log(self.init_omega)), lower=True )

        return self

    # all getters that do something
    @property
    def ell(self):
        return np.exp(self.Ll @ self.wl_ell)
    @property
    def sigma(self):
        return np.exp(self.Ls @ self.wl_sigma)
    @property
    def omega(self):
        return np.exp(self.Lo @ self.wl_omega)
    @property
    def l_ell(self):
        return self.Ll @ self.wl_ell
    @property
    def l_sigma(self):
        return self.Ls @ self.wl_sigma
    @property
    def l_omega(self):
        return self.Lo @ self.wl_omega
    # mu is always fixed to mean of that parameters' values
    @property
    def l_muell(self):
        return np.mean(self.l_ell, axis=0)
    @property
    def l_musigma(self):
        return np.mean(self.l_sigma, axis=0)
    @property
    def l_muomega(self):
        return np.mean(self.l_omega, axis=0)


#//#
#//# gpposterior.m
#//#

def gpposterior(x, y, xt, my, ell, sigma, omega):
# non-stationary (scalar) gaussian kernel

    Ktt = gausskernel(x,x,ell,sigma, omega)
    Kts = gausskernel(x,xt,ell,sigma, 0)
    Kss = gausskernel(xt,xt,ell,sigma, 0)
    Kst = Kts.T

    A = MLslash(Kst, Ktt)

    pmean = A@(y - my) + my
    pcov = Kss - A@Kts

    return pmean, pcov


#//#
#//# gradient.m
#//#

def gradient(gp, nsfuncs=None):
# Compute gradient descent over the marginal log likelihood (MLL) of the
# adaptiveGP against the 9 parameters and 3 latent functions

    if nsfuncs is None:
        nsfuncs = gp.nsfuncs

    # initial step size
    step = 1e-5

    # MLL of the initial values
    mlls = np.zeros((gp.graditers,gp.p))
    mlls[0], _, _, _, _ = nsgpmll(gp)

    if gp.verbose:
        print('%5sgp %6.d %8.2f %9.2f' % (nsfuncs, 1, np.log10(step), np.mean(mlls[0])))

    if gp.plotiters:
        plotnsgp(gp,1,1)
        MLdrawnow()

#    dml = dms = dmo = 0
    dbl = dbs = dbo = 0
    dal = das = dao = 0

    # gradient steps over all parameters
    for iter in range(1, gp.graditers):
        # derivatives of latent functions and means
        dwl_l =   deriv_ell(gp, 'l' not in nsfuncs)
        dwl_s = deriv_sigma(gp, 's' not in nsfuncs)
        dwl_o = deriv_omega(gp, 'o' not in nsfuncs)

        # deriv means
#        if 'm' in pars.nsfuncs:
#            dml, dms, dmo = deriv_mus(pars)
        if 'b' in gp.nsfuncs:
            dbl, dbs, dbo = deriv_betas(gp)
        if 'a' in gp.nsfuncs:
            dal, das, dao = deriv_alphas(gp)

        # save old parameters
        l_cp = gp.wl_ell
        s_cp = gp.wl_sigma
        o_cp = gp.wl_omega

        # gradient steps
        gp.wl_ell   = gp.wl_ell   + step*dwl_l  #//# don't use `+=` because that would modify the instance referenced by `l_cp`, which we want to preserve
        gp.wl_omega = gp.wl_omega + step*dwl_o
        gp.wl_sigma = gp.wl_sigma + step*dwl_s

        # gradients of alpha
        gp.alphaell   += step*dal
        gp.alphasigma += step*das
        gp.alphaomega += step*dao

        # gradients of beta
        gp.betaell   += step*dbl
        gp.betasigma += step*dbs
        gp.betaomega += step*dbo

        # compute MLL
        mlls[iter], _, _, _, _ = nsgpmll(gp)

        # update step
        if np.all(mlls[iter] < mlls[iter-1]):   # if overshooting, go back and decrease step size
            gp.wl_ell = l_cp
            gp.wl_sigma = s_cp
            gp.wl_omega = o_cp
            mlls[iter] = mlls[iter-1]

            step *= 0.70 # drop 25% if failing
        else:
            step *= 1.10 # increase 5% if going nicely

        if gp.verbose and ((1+iter) % 100) == 0:
            print('%5sgp %6.d %8.2f %9.2f' % (nsfuncs, 1+iter, np.log10(step), np.mean(mlls[iter])))

        if (np.log10(step) < -7) or (iter > 50 and (np.mean(mlls[iter]-mlls[iter-30]) < 0.1)):
            print('%5sgp %6.d %8.2f %9.2f' % (nsfuncs, 1+iter, np.log10(step), np.mean(mlls[iter])))
            break

        if gp.plotiters and (iter % 10) == 0:
            plotnsgp(gp,1,1)
            MLdrawnow()

    return gp


#//#
#//# logmvnpdf.m
#//#

def logmvnpdf(x,m,S):
# outputs log likelihood array for observations x  where x_n ~ N(mu,Sigma)
# x is NxD, mu is 1xD, Sigma is DxD

    logdet = lambda A: (2*np.sum(np.log(np.diag(cholesky(A, lower=False)))))

    D = MLlength(x)
    const = -0.5 * D * np.log(2*np.pi)

    term1 = -0.5 * MLslash((x-m).T, S) @ (x-m)

    try:
        term2 = const - 0.5 * logdet(S)
    except(ValueError):
        print('chol failed: increasing diagonal to fix')
        fix = np.abs(np.min(eigvals(S)))
#        np.min(np.eigvals( S + fix  ))
        term2 = const - 0.5 * logdet(S + (fix*1.01)*np.eye(*S.shape))

    return term1 + term2


#//#
#//# makesample.m
#//#

def makesample(gp):

    theta = np.empty(0)
    if 'l' in gp.nsfuncs:
        theta = np.append(theta, gp.wl_ell)
    else:
        theta = np.append(theta, gp.wl_ell[0])
    if 's' in gp.nsfuncs:
        theta = np.append(theta, gp.wl_sigma)
    else:
        theta = np.append(theta, gp.wl_sigma[0])
    if 'o' in gp.nsfuncs:
        theta = np.append(theta, gp.wl_omega)
    else:
        theta = np.append(theta, gp.wl_omega[0])

    if 'm' in gp.nsfuncs:
        theta = np.append(theta, (gp.l_muell, gp.l_musigma, gp.l_muomega))
    if 'a' in gp.nsfuncs: theta = np.append(theta, (gp.alphaell, gp.alphasigma, gp.alphaomega))
    if 'b' in gp.nsfuncs: theta = np.append(theta, (gp.betaell, gp.betasigma, gp.betaomega))

    return theta


#//#
#//# nsgausskernel.m
#//#

def nsgausskernel(x1, x2, l1, l2, s1, s2, o):
# non-stationary (scalar) gaussian kernel

    l1 = np.exp(l1)
    l2 = np.exp(l2)
    s1 = np.exp(s1)
    s2 = np.exp(s2)
    o = np.exp(o)

    sumcross = lambda v1, v2: v1[:,np.newaxis] + v2[np.newaxis]

    # loop through all inputs if variables are matrices
    K = np.outer(s1,s2) * np.sqrt((2*np.outer(l1,l2))/sumcross(l1**2,l2**2)) * np.exp(- (cdist(x1,x2,'sqeuclidean') / sumcross(l1**2,l2**2)))
    return K + np.diagflat(o**2)


#//#
#//# nsgausskernelderiv.m
#//#

def nsgausskernelderiv(x1, x2, l1, l2, s1, s2):
# ns-gauss kernel's first derivative

    sumcross = lambda v1, v2: v1[:,np.newaxis] + v2[np.newaxis]
    Kf = nsgausskernel(x1,x2,l1,l2,s1,s2,-np.inf)

    D = x1[:,np.newaxis] - x2[np.newaxis]
    L = sumcross(np.exp(l1)**2,np.exp(l2)**2)

    return - 2 * D / L * Kf


#//#
#//# nsgausskernelderiv2.m
#//#

def nsgausskernelderiv2(x1, x2, l1, l2, s1, s2):
# ns-gauss kernel's second derivative

    sumcross = lambda v1, v2: v1[:,np.newaxis] + v2[np.newaxis]
    Kf = nsgausskernel(x1,x2,l1,l2,s1,s2,-np.inf)

    D = cdist(x1,x2,'sqeuclidean')
    L = sumcross(np.exp(l1)**2,np.exp(l2)**2)

    return (2/L - 4*(D*L**(-2))) * Kf


#//#
#//# nsgp.m
#//#

def nsgp(x, y, nsfuncs='lso', optim='grad', gp=None, **kwargs):
# Main function to learn the adaptive GP
#
# INPUTS
#  - x       : input data vector
#  - y       : output data vector
#  - nsfuncs : string of components set as nonstationary [default 'lso']
#      l : lengthscale
#      s : signal variance
#      o : noise variance
#  - optim : either 'grad' [default] or 'hmc'
#  - gp : instance of gpmodel class, or None
#  - kwargs : define parameters
#
# OUTPUTS
#  - params  : parameter values learned [grad] or initialised [hmc]
#  - samples : hmc samples of the posterior? [hmc] or empty [grad]
#

    # no params given
    if not gp:
        gp = gpmodel(x,y)

    # set all remaining arguments
    gp.__dict__.update(kwargs)

    gp.nsfuncs = nsfuncs
    gp.optim = optim
    gp.init()

    samples = np.empty(0)
    if optim == 'hmc':
        # HMC sampling, do 1 chain (call this function repeatedly for more chains)
        samples,mll = nsgpnuts(gp)

    elif optim == 'grad':
        # gradient descent (several restarts included)
        gp,mll = nsgpgrad(gp)

#    mse = np.nan
#    nmse = np.nan
#    nlpd = np.nan
#    if cv < 1:
#        mse,nmse,nlpd = testerrors(gp, samples)
#        print('mse=%.4f nmse=%.4f mlpd=%.4f' % (mse, nmse, nlpd))

    mse = nmse = nlpd = None
    return gp, samples, mll, mse, nmse, nlpd


#//#
#//# nsgpgrad.m
#//#

def nsgpgrad(gp):
# learn a nonstationary GP

    print('Optimizing for %d restarts ...' % gp.restarts)
    print('  model   iter stepsize       mll')

    # perform gradient search using 'gp.restarts' initial conditions
    # the first initial condition is always fixed to the 'gp.init_ell', etc
    # remaining are randomised

    # store models
    gps = [None]*gp.restarts
    mlls = np.zeros((gp.restarts,1))

    gps[0] = gradient(gp)
    mlls[0], _, _, _, _ = nsgpmll(gp)

    for iter in range(1,gp.restarts):
        # set random initial values
        gp.init_ell = random.uniform(0.03, 0.3)
        gp.init_sigma = random.uniform(0.1, 0.5)
        gp.init_omega = random.uniform(0.01, 0.10)
        gp.init()  # do choleskies, etc.

        gps[iter] = gradient(gp)
        mlls[iter], _, _, _, _ = nsgpmll(gp)

    mll, i = np.max(mlls), np.argmax(mlls)
    gp = gps[i]

    print('Best model mll=%.2f' % mll)

#
#    gp = gradient(gp)
#    mll = nsgpmll(gp)
#    return gp, mll
#
#    bestgp = []
#    bestmll = -np.inf
#
#    if gp.nsfuncs == 'lso':
#
#        # directly go to full model
#        gp = gradient(gp)
#        mll = nsgpmll(gp)
#        if mll > bestmll:
#            bestmll = mll
#            bestgp = gp
#
#        # iteratively add more stuff
#        ords = list(reversed(list(permutations(gp.nsfuncs))))
#        for j in range(len(ords)):
#            gp = gradient(gp, '') # stationary
#            for i in range(len(gp.nsfuncs)):
#                gp = gradient(gp, ''.join(ords[j][:i+1]))
#            mll = nsgpmll(gp)
#
#            if mll > bestmll
#                bestmll = mll
#                bestgp = gp
#
#        # go from stationary to full
#        gp = gradient(gp, '')
#        gp = gradient(gp)
#        mll = nsgpmll(gp)
#        if mll > bestmll:
#            bestmll = mll
#            bestgp = gp
#
#        gp = bestgp
#        mll = bestmll
#
#        return gp, mll
#
#
#
#    ## find best standard GP solution from multiple restarts
#    # generate initial values
#    ls = np.log(np.append(0.10, random.uniform(0.02, 0.20, gp.restarts-1)))
#    ss = np.log(np.append(0.5*np.max(gp.ytr, axis=0), random.uniform(0.20, 0.80, gp.restarts-1)))
#    os = np.log(np.append(np.mean(np.abs(np.diff(gp.ytr, axis=0)), axis=0), random.uniform(0.01, 0.20, gp.restarts-1)))
#
#    n = gp.n
#
#    bestmll = -np.inf
#    bestpars = gp
#    for r in range(gp.restarts):
# #        print(' [%d/%d] Gradients of stationary parameters' % (1+r, pars.restarts))
#
# #        pars.muell = ls[r]
# #        pars.musigma = ss[r]
# #        pars.muomega = os[r]
#
#          gp.ell = ls[r]*np.ones((n,1))
#          gp.sigma = ss[r]*np.ones((n,1))
#          gp.omega = os[r]*np.ones((n,1))
#
#          if MLisfield(gp, 'Ll'):
#              gp.ell = MLbackslash(gp.Ll,(ls[r]*np.ones((n,1))))
#          if MLisfield(gp, 'Ls'):
#              gp.sigma = MLbackslash(gp.Ls,(ss[r]*np.ones((n,1))))
#          if MLisfield(gp, 'Lo'):
#              gp.omega = MLbackslash(gp.Lo,(os[r]*np.ones((n,1))))
#
#          gp = gradient(gp, '')
#
#          mll = nsgpmll(gp)
#
#          if mll > bestmll:
#              bestmll = mll
#              bestpars = gp
#
#    gp = bestpars
#    pars1 = gp
#    pars2 = gp
#
#    ## (1) add nonstationary functions one-by-one
#    for k in range(1, pars1.maxrounds+1):
#        # optimise functions one-by-one in order (several rounds)
#        for i in range(len(pars1.nsfuncs)):
# #            print(' [%d/%d] Gradients of non-stationary function "%s"' % (k, pars1.maxrounds, pars1.nsfuncs[i]))
#            pars1 = gradient(pars1, pars1.nsfuncs[i])
#
#        # always do stationary learning at end
# #        print(' Gradients of stationary function')
#        pars1 = gradient(pars1, '')
#
#        # full learning at end
# #        print(' [1/1] Gradients of non-stationary function "%s"' % pars1.nsfuncs)
#        pars1 = gradient(pars1)
#        mll1 = nsgpmll(pars1)
#
#        ## (2) do all nonstationary functions simultaneously
#        for k in range(1, pars2.maxrounds+1):
# #            print(' [%d/%d] Gradients of non-stationary function "%s"' % (k, pars2.maxrounds, pars2.nsfuncs))
#            pars2 = gradient(pars2)
#
#            # always do stationary learning at end
# #            print(' Gradients of stationary function')
#            pars2 = gradient(pars2, '')
#        mll2 = nsgpmll(pars2)
#
#        ## choose better
#        if mll1 > mll2:
#            gp = pars1
#            mll = mll1
#        else:
#            gp = pars2
#            mll = mll2
#
#    print('Best model mll=%.2f' % mll)
#

    return gp, mll


#//#
#//# nsgpmll.m
#//#

def nsgpmll(pars):

    n = MLlength(pars.xtr)
    if ('a' in pars.nsfuncs) or ('b' in pars.nsfuncs):
        pars.Kl = gausskernel(pars.xtr, pars.xtr, pars.betaell,   pars.alphaell,   pars.tol)
        pars.Ks = gausskernel(pars.xtr, pars.xtr, pars.betasigma, pars.alphasigma, pars.tol)
        pars.Ko = gausskernel(pars.xtr, pars.xtr, pars.betaomega, pars.alphaomega, pars.tol)
    Ky = nsgausskernel(pars.xtr, pars.xtr, pars.l_ell, pars.l_ell, pars.l_sigma, pars.l_sigma, pars.l_omega)
    zs = np.zeros((pars.n,pars.p))

    # check if non-sdp or low condition number
    try:
        cholesky(Ky, lower=False)
        p = 0
    except(ValueError):
        p = 1
    rc = 1. / cond(Ky, p=1)
    if p > 0 or rc < 1e-15:
        val = -np.inf
        return val, None, None, None, None

    # assuming exp-transformation here
    valy = np.diag(logmvnpdf(pars.ytr, zs, Ky))

    if MLlength(pars.ell) == 1:
        vall = logmvnpdf(pars.ell*np.ones((n,1)), muell, pars.Kl)
    else:
        vall = logmvnpdf(pars.l_ell, pars.l_muell, pars.Kl)

    if MLlength(pars.sigma) == 1:
        vals = logmvnpdf(pars.sigma*np.ones((n,1)), musigma, pars.Ks)
    else:
        vals = logmvnpdf(pars.l_sigma, pars.l_musigma, pars.Ks)

    if MLlength(pars.omega) == 1:
        vale = logmvnpdf(pars.omega*np.ones((n,1)), muomega, pars.Ko)
    else:
        vale = logmvnpdf(pars.l_omega, pars.l_muomega, pars.Ko)

    val = valy + vall + vals + vale

    return val, valy, vall, vals, vale


#//#
#//# nsgpnuts.m
#//#

def nsgpnuts(pars):
# HMC-NUTS-DA sampling with Hoffman/Gelman code
# hmc parameters:
#  pars.nsfuncs is a string of derivatives to include, defaults to 'sle'
#    where 'l' = ell   [lengthscale]
#          's' = sigma [signal]
#          'o' = omega [error]
#          'm' = prior means
#          'a' = prior alphas
#          'b' = prior betas
#  pars.warmups  : default 1000
#  pars.iters    : default 1000
#  pars.maxdepth : default 9
#  pars.delta    : default 0.5
#

    def gradients(theta):
        nonlocal iter, pars

        # put sample into pars
        pars = parsesample(theta,pars)

        # compute logp and grads
        logp, _, _, _, _ = nsgpmll(pars)
        grad = computegrads(pars)

        if pars.plotiters and iter > 20 and (iter%100)==0:
            plotnsgp(pars,1,1)
            MLdrawnow()

        # iter counter
        iter += 1

        return logp, grad

    ## defaults if not specified before
    if not MLisfield(pars, 'iters'):
        pars.iters = 1000
    if not MLisfield(pars, 'warmups'):
        pars.warmups = 1000
    if not MLisfield(pars, 'delta'):
        pars.delta = 0.5
    if not MLisfield(pars, 'nsfuncs'):
        pars.nsfuncs = 'slo'

    f = gradients
    iter = 1

    # take MAP optimal as starting point
    # parsstat = gradient(xtr, ytr, pars, '')
    # parsmap = gradient(xtr, ytr, parsstat)

    # make initial sample
    theta0 = makesample(pars)

    # compute HMC-NUTS
    if 'l' in pars.nsfuncs:
        samples = nuts(pars.leapfrog, f, pars.iters, theta0)
    else:
        samples, _ = nuts_da(f, pars.iters, pars.warmups, theta0, pars.delta, pars.maxdepth, pars.verbose)

    # compute MLLS
    mlls = np.zeros((pars.iters,1))
    for i in range(pars.iters):
        pars = parsesample(samples[i], pars)
        mlls[i], _, _, _, _ = nsgpmll(pars)

    return samples, mlls


#//#
#//# nsgpposterior.m
#//#

def nsgpposterior(gp, xt=None): #, computederiv=False):
# Model posteriors over target points 'xt'
# returns:
#  - fmean      : f posterior mean
#  - fstd       : f posterior std
#  - lt         :   ell posterior mean [lengthscale]
#  - st         : sigma posterior mean [signal variance]
#  - ot         : omega posterior mean [noise variance]
#  - fderivmean : f derivative posterior mean
#  - fderivstd  : f derivative posterior std

    if xt is None:
        xt = gp.xtr

    ot, _ = gpposterior(gp.xtr, gp.l_omega, xt, gp.l_muomega, gp.betaomega, gp.alphaomega, gp.tol)
    lt, _ = gpposterior(gp.xtr, gp.l_ell,   xt, gp.l_muell,   gp.betaell,   gp.alphaell,   gp.tol)
    st, _ = gpposterior(gp.xtr, gp.l_sigma, xt, gp.l_musigma, gp.betasigma, gp.alphasigma, gp.tol)

    # extrapolate unknown function
    Ktt = nsgausskernel(gp.xtr, gp.xtr, gp.l_ell, gp.l_ell,   gp.l_sigma, gp.l_sigma, gp.l_omega)
    Kts = nsgausskernel(gp.xtr, xt,     gp.l_ell, lt,         gp.l_sigma, st,         -np.inf)
    Kss = nsgausskernel(xt,     xt,     lt,       lt,         st,         st,         -np.inf)
    Kst = Kts.T

    A = MLslash(Kst, Ktt)

    fmean = A@(gp.ytr - gp.muf) + gp.muf
    fcov = Kss - A@Kts
    fstd = np.sqrt(np.diag(fcov))

    # compute also the derivative GP
#    if computederiv:
#        Kdst = nsgausskernelderiv(xt,  gp.xtr, lt, gp.l_ell, st, gp.l_sigma)
#        Kdss = nsgausskernelderiv2(xt, xt,      lt, lt,  st, st)
#        fderivmean = MLslash(Kdst, Ktt) @ gp.ytr
#        fderivstd = np.sqrt(np.diag(Kdss - MLslash(Kdst, Ktt) @ Kdst.T))
#
#        _,_,fmean,fstd,lt,st,ot,fderivmean,fderivstd = denormalise(gp, fmean, fstd, np.exp(lt), np.exp(st), np.exp(ot), fderivmean, fderivstd)

    lt = np.exp(lt)
    st = np.exp(st)
    ot = np.exp(ot)

    _,_,fmean,fstd,lt,st,ot = denormalise(gp, fmean, fstd, lt, st, ot)

#    return fmean, fstd, lt, st, ot, fderivmean, fderivstd
    return fmean, fstd, lt, st, ot


#//#
#//# nuts.m
#//#

def nuts(epsilon, f, M, theta0):
# samples = nuts(epsilon, f, M, theta0)
#
# Implements the No-U-Turn Sampler (NUTS), specifically, algorithm 3
# from the NUTS paper (Hoffman & Gelman, 2011).
#
# This algorithm requires that you specify a step size parameter epsilon---try
# nuts_da() if you want to avoid tuning epsilon.
#
# epsilon is a step size parameter.
# f(theta) should be a function that returns the log probability its
# gradient evaluated at theta. I.e., you should be able to call
# logp, grad = f(theta).
# M is the number of samples to generate.
# theta0 is a 1-by-D vector with the desired initial setting of the parameters.
#
# The returned variable "samples" is an M-by-D matrix of samples generated
# by NUTS. samples(1, :) = theta0.

# Copyright (c) 2011-2012, Matthew D. Hoffman
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    def leapfrog(theta, r, grad, epsilon, f):
        rprime = r + 0.5 * epsilon * grad
        thetaprime = theta + epsilon * rprime
        logpprime, gradprime = f(thetaprime)
        rprime += 0.5 * epsilon * gradprime
        return thetaprime, rprime, gradprime, logpprime

    def stop_criterion(thetaminus, thetaplus, rminus, rplus):
        thetavec = thetaplus - thetaminus
        return (thetavec @ rminus >= 0) and (thetavec @ rplus >= 0)

    # The main recursion.
    def build_tree(theta, r, grad, logu, v, j, epsilon, f):
        if (j == 0):
            # Base case: Take a single leapfrog step in the direction v.
            thetaprime, rprime, gradprime, logpprime = leapfrog(theta, r, grad, v*epsilon, f)
            joint = logpprime - 0.5 * (rprime @ rprime)
            # Is the new point in the slice?
            nprime = np.array(logu < joint, dtype=int)
            # Is the simulation wildly inaccurate?
            sprime = logu - 1000 < joint
            # Set the return values---minus=plus for all things here, since the
            # "tree" is of depth 0.
            thetaminus = thetaprime
            thetaplus = thetaprime
            rminus = rprime
            rplus = rprime
            gradminus = gradprime
            gradplus = gradprime
        else:
            # Recursion: Implicitly build the height j-1 left and right subtrees.
            thetaminus, rminus, gradminus, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime = \
                build_tree(theta, r, grad, logu, v, j-1, epsilon, f)
            # No need to keep going if the stopping criteria were met in the first
            # subtree.
            if sprime:
                if (v == -1):
                    thetaminus, rminus, gradminus, _, _, _, thetaprime2, gradprime2, logpprime2, nprime2, sprime2 = \
                        build_tree(thetaminus, rminus, gradminus, logu, v, j-1, epsilon, f)
                else:
                    _, _, _, thetaplus, rplus, gradplus, thetaprime2, gradprime2, logpprime2, nprime2, sprime2 = \
                        build_tree(thetaplus, rplus, gradplus, logu, v, j-1, epsilon, f)
                # Choose which subtree to propagate a sample up from.
                if (random.uniform() < nprime2 / (nprime + nprime2 + 1.e-40)):  #//# in MATLAB division by zero gives NaN, and :mat:`(anything < NaN)` evaluates to 0; we need to mimic this without actually dividing by zero
                    thetaprime = thetaprime2
                    gradprime = gradprime2
                    logpprime = logpprime2
                # Update the number of valid points.
                nprime += nprime2
                # Update the stopping criterion.
                sprime &= sprime2 and stop_criterion(thetaminus, thetaplus, rminus, rplus)
        return thetaminus, rminus, gradminus, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime

    assert theta0.ndim == 1

    D = MLlength(theta0)
    samples = np.zeros((M, D))

    logp,grad = f(theta0)
    samples[0] = theta0

    for m in range(1,M):
        #//#print('sample %4.d  logp %.3f' % (1+m, logp))
        print('sample', '%4.d' % (1+m), ' logp', *( '%.3f' % _ for _ in logp.ravel() ))

        # resample momenta
        r0 = random.normal(size=D)
        # Joint log-probability of theta and momenta r.
        joint = logp - 0.5 * (r0 @ r0)
        # Resample u ~ uniform([0, exp(joint)]).
        # Equivalent to (log(u) - joint) ~ exponential(1).
        logu = joint - random.exponential(1)
        # Initialize tree.
        thetaminus = samples[m-1]
        thetaplus = samples[m-1]
        rminus = r0
        rplus = r0
        gradminus = grad
        gradplus = grad
        # Initial height j = 0.
        j = 0
        # If all else fails, the next sample is the previous sample.
        samples[m] = samples[m-1]
        # Initially the only valid point is the initial point.
        n = 1

        # Main loop---keep going until the criterion s == 0.
        s = True
        while s:
            # Choose a direction. -1=backwards, 1=forwards.
            v = 2*(random.uniform() < 0.5)-1
            # Double the size of the tree.
            if (v == -1):
                thetaminus, rminus, gradminus, _, _, _, thetaprime, gradprime, logpprime, nprime, sprime = \
                    build_tree(thetaminus, rminus, gradminus, logu, v, j, epsilon, f)
            else:
                _, _, _, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime = \
                    build_tree(thetaplus, rplus, gradplus, logu, v, j, epsilon, f)
            # Use Metropolis-Hastings to decide whether or not to move to a
            # point from the half-tree we just generated.
            if (sprime and (random.uniform() < nprime/n)):
                samples[m] = thetaprime
                logp = logpprime
                grad = gradprime
            # Update number of valid points we've seen.
            n += nprime
            # Decide if it's time to stop.
            s = sprime and stop_criterion(thetaminus, thetaplus, rminus, rplus) and j < 9
            # Increment depth.
            j += 1

    return samples


#//#
#//# nuts_da.m
#//#

def nuts_da(f, M, Madapt, theta0, delta=0.6, maxdepth=10, verbose=False):
# samples, epsilon = nuts_da(f, M, Madapt, theta0, delta)
#
# Implements the No-U-Turn Sampler (NUTS), specifically, algorithm 6
# from the NUTS paper (Hoffman & Gelman, 2011). Runs Madapt steps of
# burn-in, during which it adapts the step size parameter epsilon, then
# starts generating samples to return.
#
# epsilon is a step size parameter.
# f(theta) should be a function that returns the log probability its
# gradient evaluated at theta. I.e., you should be able to call
# logp, grad = f(theta).
# M is the number of samples to generate.
# Madapt is the number of steps of burn-in/how long to run the dual
# averaging algorithm to fit the step size epsilon.
# theta0 is a 1-by-D vector with the desired initial setting of the parameters.
# delta should be between 0 and 1, and is a target HMC acceptance
# probability. Defaults to 0.6 if unspecified.
#
# The returned variable "samples" is an M-by-D matrix of post-burn-in samples
# generated by NUTS.
# The returned variable "epsilon" is the step size that was fit using dual
# averaging.

# Copyright (c) 2011, Matthew D. Hoffman
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    global nfevals
    nfevals = 0

    def leapfrog(theta, r, grad, epsilon, f):
        rprime = r + 0.5 * epsilon * grad
        thetaprime = theta + epsilon * rprime
        logpprime, gradprime = f(thetaprime)
        rprime += 0.5 * epsilon * gradprime
        global nfevals
        nfevals += 1
        return thetaprime, rprime, gradprime, logpprime

    def stop_criterion(thetaminus, thetaplus, rminus, rplus):
        thetavec = thetaplus - thetaminus
        return (thetavec @ rminus >= 0) and (thetavec @ rplus >= 0)

    # The main recursion.
    def build_tree(theta, r, grad, logu, v, j, epsilon, f, joint0):
        if (j == 0):
            # Base case: Take a single leapfrog step in the direction v.
            thetaprime, rprime, gradprime, logpprime = leapfrog(theta, r, grad, v*epsilon, f)
            joint = logpprime - 0.5 * (rprime @ rprime)
            # Is the new point in the slice?
            nprime = logu < joint
            # Is the simulation wildly inaccurate?
            sprime = logu - 1000 < joint
            # Set the return values---minus=plus for all things here, since the
            # "tree" is of depth 0.
            thetaminus = thetaprime
            thetaplus = thetaprime
            rminus = rprime
            rplus = rprime
            gradminus = gradprime
            gradplus = gradprime
            # Compute the acceptance probability.
            alphaprime = min(1, np.exp(logpprime - 0.5 * (rprime @ rprime) - joint0))
            nalphaprime = 1
        else:
            # Recursion: Implicitly build the height j-1 left and right subtrees.
            thetaminus, rminus, gradminus, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime, alphaprime, nalphaprime = \
                build_tree(theta, r, grad, logu, v, j-1, epsilon, f, joint0)
            # No need to keep going if the stopping criteria were met in the first
            # subtree.
            if sprime:
                if (v == -1):
                    thetaminus, rminus, gradminus, _, _, _, thetaprime2, gradprime2, logpprime2, nprime2, sprime2, alphaprime2, nalphaprime2 = \
                        build_tree(thetaminus, rminus, gradminus, logu, v, j-1, epsilon, f, joint0)
                else:
                    _, _, _, thetaplus, rplus, gradplus, thetaprime2, gradprime2, logpprime2, nprime2, sprime2, alphaprime2, nalphaprime2 = \
                        build_tree(thetaplus, rplus, gradplus, logu, v, j-1, epsilon, f, joint0)
                # Choose which subtree to propagate a sample up from.
                if (random.uniform() < nprime2 / (nprime + nprime2 + 1.e-40)):  #//# in MATLAB division by zero gives NaN, and :mat:`(anything < NaN)` evaluates to 0; we need to mimic this without actually dividing by zero
                    thetaprime = thetaprime2
                    gradprime = gradprime2
                    logpprime = logpprime2
                # Update the number of valid points.
                nprime += nprime2
                # Update the stopping criterion.
                sprime &= sprime2 and stop_criterion(thetaminus, thetaplus, rminus, rplus)
                # Update the acceptance probability statistics.
                alphaprime += alphaprime2
                nalphaprime += nalphaprime2
        return thetaminus, rminus, gradminus, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime, alphaprime, nalphaprime

    def find_reasonable_epsilon(theta0, grad0, logp0, f):
        epsilon = 0.1
        r0 = random.normal(size=MLlength(theta0))
        # Figure out what direction we should be moving epsilon.
        _, rprime, _, logpprime = leapfrog(theta0, r0, grad0, epsilon, f)
        acceptprob = np.exp(logpprime - logp0 - 0.5 * (rprime @ rprime - r0 @ r0))
        a = 2 * (acceptprob > 0.5) - 1
        # Keep moving epsilon in that direction until acceptprob crosses 0.5.
        while (acceptprob**a > 2**(-a)):
            epsilon *= 2**a
            _, rprime, _, logpprime = leapfrog(theta0, r0, grad0, epsilon, f)
            acceptprob = np.exp(logpprime - logp0 - 0.5 * (rprime @ rprime - r0 @ r0))
        return epsilon

    assert theta0.ndim == 1

    D = MLlength(theta0)
    samples = np.zeros((M+Madapt, D))

    logp, grad = f(theta0)
    samples[0] = theta0

    # Choose a reasonable first epsilon by a simple heuristic.
    epsilon = find_reasonable_epsilon(theta0, grad, logp, f)

    if verbose:
        print('found reasonable epsilon %.5f' % epsilon)

    # Parameters to the dual averaging algorithm.
    gamma = 0.05
    t0 = 10
    kappa = 0.75
    mu = np.log(10*epsilon)
    # Initialize dual averaging algorithm.
    epsilonbar = 0.01
    Hbar = 0

    for m in range(1,M+Madapt):
        if verbose and (m % 1):
            print('sample=%d  log(epsilon)=%.1f' % (1+m, np.log10(epsilon)))
        # Resample momenta.
        r0 = random.normal(size=D)
        # Joint log-probability of theta and momenta r.
        joint = logp - 0.5 * (r0 @ r0)
        # Resample u ~ uniform([0, exp(joint)]).
        # Equivalent to (log(u) - joint) ~ exponential(1).
        logu = joint - random.exponential(1)
        # Initialize tree.
        thetaminus = samples[m-1]
        thetaplus = samples[m-1]
        rminus = r0
        rplus = r0
        gradminus = grad
        gradplus = grad
        # Initial height j = 0.
        j = 0
        # If all else fails, the next sample is the previous sample.
        samples[m] = samples[m-1]
        # Initially the only valid point is the initial point.
        n = 1

        # Main loop---keep going until the criterion s == 0.
        s = True
        while s:

#            if verbose:
#                print(' depth=%d logp=%.2f' % (j, logp))

            # Choose a direction. -1=backwards, 1=forwards.
            v = 2*(randon.uniform() < 0.5)-1
            # Double the size of the tree.
            if (v == -1):
                thetaminus, rminus, gradminus, _, _, _, thetaprime, gradprime, logpprime, nprime, sprime, alpha, nalpha = \
                    build_tree(thetaminus, rminus, gradminus, logu, v, j, epsilon, f, joint)
            else:
                _, _, _, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime, alpha, nalpha = \
                    build_tree(thetaplus, rplus, gradplus, logu, v, j, epsilon, f, joint)
            # Use Metropolis-Hastings to decide whether or not to move to a
            # point from the half-tree we just generated.
            if (sprime and (random.uniform() < nprime/n)):
                samples[m] = thetaprime
                logp = logpprime
                grad = gradprime
            # Update number of valid points we've seen.
            n += nprime
            # Decide if it's time to stop.
#            s = sprime and stop_criterion(thetaminus, thetaplus, rminus, rplus)
                # MARKUS: put a maximum tree depth in
            s = sprime and stop_criterion(thetaminus, thetaplus, rminus, rplus) and j < maxdepth
            # Increment depth.
            j += 1

#        print(' logp=%.2f' % logp)

        # Do adaptation of epsilon if we're still doing burn-in.
        eta = 1 / (m - 1 + t0)
        Hbar = (1 - eta) * Hbar + eta * (delta - alpha / nalpha)
        if (m <= Madapt):
            epsilon = np.exp(mu - np.sqrt(m-1)/gamma * Hbar)
            eta = (m-1)**-kappa
            epsilonbar = np.exp((1 - eta) * np.log(epsilonbar) + eta * np.log(epsilon))
        else:
            epsilon = epsilonbar

    samples = samples[Madapt:]
    print('Took %d gradient evaluations.\n' % nfevals)

    return samples, epsilon


#//#
#//# parsesample.m
#//#

def parsesample(theta, pars):

    # retrieve results from 'theta' into 'pars'

    n = MLlength(pars.ell)
    j = 0
    if 'l' in pars.nsfuncs:
        pars.wl_ell   = theta[j:j+n]
        j += n
    else:
        pars.wl_ell   = np.full(n, theta[j])
        j += 1
    if 's' in pars.nsfuncs:
        pars.wl_sigma = theta[j:j+n]
        j += n
    else:
        pars.wl_sigma = np.full(n, theta[j])
        j += 1
    if 'o' in pars.nsfuncs:
        pars.wl_omega = theta[j:j+n]
        j += n
    else:
        pars.wl_omega = np.full(n, theta[j])
        j += 1
    if 'm' in pars.nsfuncs:
        pars.muell   = theta[j]
        pars.musigma = theta[j+1]
        pars.muomega = theta[j+2]
        j += 3
    if 'a' in pars.nsfuncs:
        pars.alphaell   = theta[j]
        pars.alphasigma = theta[j+1]
        pars.alphaomega = theta[j+2]
        j += 3
    if 'b' in pars.nsfuncs:
        pars.betaell   = theta[j]
        pars.betasigma = theta[j+1]
        pars.betaomega = theta[j+2]

    return pars


#//#
#//# plotgp.m
#//#

@import_plt
def plotgp(x, y, xt, pars):

    # MLL of the initial values
    mll = gpmll(x,y,pars)

    xt = np.linspace(0, 25, 200)[:,np.newaxis]
    nt = MLlength(xt)
    yt,ytcov = gpposterior(x, y, xt, 0, pars.ell, pars.sigma, pars.eps, pars.tol)
    ytvar = np.diag(ytcov)
    ytstd = np.sqrt(ytvar)

    plt.subplot(2,1,1)
    h = plt.fill(np.append(xt, np.flip(xt)), np.append(yt[:,0]+2*ytstd, np.flip(yt[:,0]-2*ytstd)), 'red', alpha=0.20, edgecolor=None)
    plt.plot(x,y,'ko')
    plt.plot(xt,yt, 'k-', linewidth=1.5)
    plt.plot(xt, yt + 2*np.sqrt(ytstd**2 + pars.eps**2), 'k--')
    plt.plot(xt, yt - 2*np.sqrt(ytstd**2 + pars.eps**2), 'k--')
    plt.title('Stationary GP posterior [MLL  %.3f]' % mll)
    plt.xlabel('time')
    plt.ylabel('value')
    plt.legend(['95% posterior','training data','mean posterior','95% noisy posterior'])

    plt.subplot(2,1,2)
    K = gausskernel(xt,xt,pars.ell,pars.sigma,0)
    Z = random.multivariate_normal(np.zeros(nt), K,5).T
    plt.plot(xt,Z, '-x')
    plt.title('5 samples from the prior N(0,K)')


#//#
#//# plotnsgp.m
#//#

@import_plt
def plotnsgp(gp, plotlatent=False, plotderivs=False, plotkernel=False, ts=False, truemodel=None):
# Plots the GP model
# parameters:
#  - plotlatent : plot latent functions (ell,sigma,omega)
#  - plotderivs : plot derivative GP
#  - plotkernel : plot covariance matrix
#  - ts         : test data to include
#  - truemodel  : highlight true function
#

    ## 2D plot
    if gp.xtr.shape[1] == 2:
        plotnsgp2d(gp)
        return

    ##

    squares = 1 + plotkernel + plotlatent + plotderivs
    cols = np.array([[248, 118, 109], [0, 186, 56], [97, 156, 255]]) / 255 # ggplot2 colors
    xt = np.linspace(0,1,250)[:,np.newaxis]

    # MLL of the initial values
    nsmll, _, _, _, _ = nsgpmll(gp)

    ft,ftstd,lt,st,ot = nsgpposterior(gp,xt)
    xt = np.linspace(0,1,250)[:,np.newaxis] * gp.xscale + gp.xbias

    xtr,ytr,_,_,_,_,_ = denormalise(gp)
#    xtr = gp.xtr
#    ytr = gp.ytr

    # test predictions
    if ts:
        mse,_,nlpd = testerrors(gp)

    ## function posterior
    if squares > 1:
        plt.subplot(squares,1,1)

    # plot 1D curve
    for i in range(gp.p):
        h = plt.fill(np.append(xt, np.flip(xt)), np.append(ft[:,i]+2*ftstd, np.flip(ft[:,i]-2*ftstd)), 'black', alpha=0.45, edgecolor=None)

    for i in range(gp.p):
        h3 = plt.fill(np.append(xt, np.flip(xt)), np.append(ft[:,i]+2*ftstd, np.flip(ft[:,i]+2*np.sqrt(ftstd**2+ot**2))), 'black', alpha=0.2, edgecolor=None)
        plt.fill(np.append(xt, np.flip(xt)), np.append(ft[:,i]-2*ftstd, np.flip(ft[:,i]-2*np.sqrt(ftstd**2+ot**2))), 'black', alpha=0.2, edgecolor=None)

    h1 = plt.plot(xtr, ytr, 'o', color='black', markersize=6)
    plt.plot(xtr, ytr, '.', color='white', markersize=12)

    if truemodel is not None:
        h4 = plt.plot(truemodel.x, truemodel.f, 'r-', linewidth=1.0)

    if ts:
        plt.title('Nonstationary function 95%% posterior [MLL=%.3f, MSE=%.3f, NLPD=%.3f]' % (nsmll,mse,nlpd))
    else:
        if (gp.nsfuncs == ''):
            plt.title('Stationary GP function 95% posterior')
        else:
            plt.title('Nonstationary %s-GP function 95%% posterior' % str.upper(gp.nsfuncs))

    plt.xlabel('time')
    plt.ylabel('value')
#    plt.legend([h1 h2 h h3], ['Data','Posterior mean','95% posterior','95% noisy posterior'], loc='upper left')
#    if truemodel is not None:
#        plt.legend([h1 h2 h h3 h4], ['Data','Posterior mean','Posterior','Noisy posterior', 'True function'], loc='upper left')
    plt.rcParams['font.size'] = 12

    ## latent functions
    if plotlatent:
        # latent functions
        plt.subplot(squares,1,2)
        plt.plot(xt, lt, color=cols[0])
        plt.plot(xt, st, color=cols[1])
        plt.plot(xt, ot, color=cols[2])

#        plt.plot(xtr, np.exp(ell),   '.', color=cols[0])
#        plt.plot(xtr, np.exp(sigma), '.', color=cols[1])
#        plt.plot(xtr, np.exp(omega), '.', color=cols[2])

        plt.xlabel('time')
        plt.ylabel('value')
        plt.ylim(0,np.max([np.max(st),np.max(lt),np.max(ot)])*1.1)
        l = plt.legend([r'$\ell$ lengthscale',r'$\sigma^2$ signal variance',r'$\omega^2$ noise variance'], loc='upper left')
        plt.rcParams['text.usetex'] = True

#        plt.legend([ r'lengthscale (ell)  [\mu=%.2f \alpha=%.2f \beta=%.2f]'       % (np.exp(gp.muell),   gp.alphaell,   gp.betaell),
#                     r'signal variance (sigma)  [\mu=%.2f \alpha=%.2f \beta=%.2f]' % (np.exp(gp.musigma), gp.alphasigma, gp.betasigma),
#                     r'noise std (omega)  [\mu=%.2f \alpha=%.2f \beta=%.2f]'       % (np.exp(gp.muomega), gp.alphaomega, gp.betaomega) ])
#        if truemodel is not None:
#            plt.legend(['Lengthscale','Signal variance','Noise variance','Generating lengthscale','Generating signal variance','Generating noise variance'])
        plt.title('Parameters')

    ## derivatives
    if plotderivs:
        plt.subplot(squares,1,3)

        dwl_l = deriv_ell(gp,   'l' not in gp.nsfuncs)
        dwl_s = deriv_sigma(gp, 's' not in gp.nsfuncs)
        dwl_o = deriv_omega(gp, 'o' not in gp.nsfuncs)

        dl_l = gp.Ll@dwl_l
        dl_s = gp.Ls@dwl_s
        dl_o = gp.Lo@dwl_o

        plt.stem(xtr, dl_o, color=cols[2])
        plt.stem(xtr, dl_s, color=cols[1])
        plt.stem(xtr, dl_l, color=cols[0])

        plt.title('Latent derivatives')

    ## kernel
    if plotkernel:
        plt.subplot(squares,1,squares)
        K = nsgausskernel(xt,xt,log(lt/gp.yscale),log(lt/gp.yscale),log(st/gp.yscale),log(st/gp.yscale),log(ot/gp.yscale))
        plt.matshow(np.real(K))
        plt.title('Kernel matrix')


#//#
#//# plotnsgp2d.m
#//#

@import_plt
def plotnsgp2d(gp):
    # grid
    ng = 70
    nl = 10
    X1,X2 = np.meshgrid(np.linspace(0,1,ng), np.linspace(0,1,ng))
    Xt = np.hstack((X1.reshape(-1,1,order='F'), X2.reshape(-1,1,order='F')))

#    # grid
#    mg = 80
#    xrng = np.max(gp.xtr, axis=0) - np.min(gp.xtr, axis=0)
#    x1 = np.linspace(np.min(gp.xtr[:,0])-0.05, np.max(gp.xtr[:,0])+0.05, mg)[:,np.newaxis]
#    x2 = np.linspace(np.min(gp.xtr[:,1])-0.05, np.max(gp.xtr[:,1])+0.05, mg)[:,np.newaxis]
#    X1,X2 = np.meshgrid(x1,x2)
#    Xt = np.hstack((X1.reshape(-1,1,order='F'), X2.reshape(-1,1,order='F')))

#    xtr,ytr,ft,ftstd,lt,st,ot,ftderiv,ftderivstd = denormalise(gp, ft, ftstd, lt, st, ot)
#    ell,sigma,omega = latentchols(gp)
#    ell = np.exp(ell)
#    sigma = np.exp(sigma)
#    omega = np.exp(omega)
#    xtr = gp.xtr
#    ytr = gp.ytr

#    ell = np.exp(ell) * gmean(gp.xscale)
#    sigma = np.exp(sigma) * gp.yscale
#    omega = np.exp(omega) * gp.yscale
#    xtr,ytr,_,_,_,_,_,_,_ = denormalise(gp)


    ft,ftstd,lt,st,ot = nsgpposterior(gp, Xt)

    # model as grid
#    ot = np.exp(gpposterior(pars.xtr, omega, Xt, pars.muomega, pars.betaomega, pars.alphaomega, pars.tol)[0])
#    lt = np.exp(gpposterior(pars.xtr, ell,   Xt, pars.muell,   pars.betaell,   pars.alphaell,   pars.tol)[0])
#    st = np.exp(gpposterior(pars.xtr, sigma, Xt, pars.musigma, pars.betasigma, pars.alphasigma, pars.tol)[0])

    plt.subplot(2,3,1)
    plt.contour(X1,X2, np.reshape(yt,(ng,ng)),nl)
    plt.title('E[f]')
    plt.colorbar()

    plt.scatter(gp.xtr[:,0], gp.xtr[:,1], 20, gp.ytr)

    plt.subplot(2,3,2)
    plt.contour(X1,X2, np.reshape(np.diag(ytcov),(ng,ng)),nl)
    plt.title('Var[f]')
    plt.colorbar()

    plt.subplot(2,3,4)
    plt.contour(X1,X2, np.reshape(lt,(ng,ng)),nl)
    plt.title('ell')
    plt.colorbar()

    plt.subplot(2,3,5)
    plt.contour(X1,X2, np.reshape(st,(ng,ng)),nl)
    plt.title('sigma')
    plt.colorbar()

    plt.subplot(2,3,6)
    plt.contour(X1,X2, np.reshape(ot,(ng,ng)),nl)
    plt.title('omega')
    plt.colorbar()


#//#
#//# plotnsgpsamples.m
#//#

@import_plt
def plotnsgpsamples(pars, samples, plotlatent=False, n=100, mapgp=None, truemodel=None):
    squares = 2 + plotlatent*3
    cols = np.array([[248, 118, 109], [0, 186, 56], [97, 156, 255]]) / 255

    # MLL of the initial values
    xt = np.linspace(0,1,200)[:,np.newaxis]
    nt = MLlength(xt)
    m = samples.shape[0]

    ## first check mll's
    mlls = np.zeros((m,1))  #//# this was `n`, but I think it should be `m` because `n < m` is anticipated below
    for i in range(m):
        pars = parsesample(samples[i], pars)
        mlls[i], _, _, _, _ = nsgpmll(pars)

    # threshold worst 10% away
    thr = np.quantile(mlls[:,0],0.10)
    samples = samples[mlls[:,0]>thr]
    m = samples.shape[0]

    ## choose indices from remaining
    if n < m:
        I = random.permutation(m)[:n]
    else:
        I = np.arange(m)

    ## precompute all functions and latents over all samples
    ots = np.zeros((n,nt))
    lts = np.zeros((n,nt))
    sts = np.zeros((n,nt))
    stds = np.zeros((n,nt))
    means = np.zeros((n,nt))
    for i in range(n):
        pars = parsesample(samples[I[i]], pars)
        for arr, _ in zip((means,stds,lts,sts,ots), nsgpposterior(pars, xt)):
            arr[i] = _.ravel()

#    xtr = pars.xtr
#    ytr = pars.ytr

    ## map GP
    if mapgp:
        ftmap,_,ltmap,stmap,otmap = nsgpposterior(mapgp, xt)
#        ftmap = ftmap * mapgp.yscale + mapgp.ybias
#        ltmap *= mapgp.yscale
#        stmap *= mapgp.yscale

    # scale
#    xt = xt * pars.xscale + pars.xbias
#    means = means * pars.yscale + pars.ybias
#    stds *= pars.yscale
#    ots *= pars.yscale
#    sts *= pars.yscale
#    xtr = pars.xtr * pars.xscale + pars.xbias
#    ytr = pars.ytr * pars.yscale + pars.ybias

#    ft,ftstd,lt,st,ot,ftderiv,ftstdderiv = nsgpposterior(pars,xt)
    xtr,ytr,_,_,_,_,_ = denormalise(pars)
#    plt.fill(np.append(xt, np.flip(xt)), np.append(ft+2*ftstd, np.flip(ft-2*ftstd)), 'black', alpha=0.25)

    ## plot function
    falpha = 0.65

    if squares > 1:
        plt.subplot(squares, 1, (1, 2))

    # plot posterior samples
    l = [None]*2
    for i in range(n):
        yt1 = means[i]+2*np.sqrt(stds[i]**2 + ots[i]**2)
        yt2 = means[i]-2*np.sqrt(stds[i]**2 + ots[i]**2)
        l[0] = plt.fill(np.append(xt, np.flip(xt)), np.append(yt1, np.flip(yt2)), 'black', alpha=falpha/n, edgecolor=None)

    if mapgp is not None:
        l.append( plt.plot(xt, ftmap, 'k--', linewidth=1.5) )
    if truemodel is not None:
        l.append( plt.plot(truemodel.x, truemodel.f, 'r-') )

    l[1] = plt.plot(xtr, ytr, 'o', color='black', markersize=6)
    plt.plot(xtr, ytr, '.', color='white', markersize=12)

    if (mapgp is not None) and (truemodel is not None):
        plt.legend(l, ['Posterior samples','Data','MAP posterior','True function'], loc='upper left')
    elif mapgp is not None:
        plt.legend(l, ['Posterior samples','Data','MAP posterior'], loc='upper left')
    elif truemodel is not None:
        plt.legend(l, ['Posterior samples','Data','True function'], loc='upper left')
    else:
        plt.legend(l, ['Posterior samples','Data'], loc='upper left')

#    plt.ylim( np.quantile(np.min(means - 2*sqrt(stds.^2 + ots.^2),axis=1), 0.05), quantile(np.max(means + 2*sqrt(stds.^2 + ots.^2),axis=1), 0.95) )
    plt.ylabel('value')

    if (pars.nsfuncs == ''):
        plt.title('Stationary GP function 95% posterior')
    else:
        plt.title('Nonstationary %s-GP function 95%% posterior' % str.upper(pars.nsfuncs))

    ## plot latent
    lw_s = 0.020
    lw_l = 0.010
    lw_o = 0.005
    alpha_s = 5
    alpha_l = 4
    alpha_o = 2

    truemodel = None

    if plotlatent:
        plt.subplot(squares,1,3)
        h = [None]
        for i in range(n):
            h[0] = plt.fill(np.append(xt, np.flip(xt)), np.append(lts[i] + lw_l, np.flip(lts[i] - lw_l)), color=cols[0], alpha=alpha_l/n, edgecolor=None)
        if mapgp is not None:
            h.append( plt.plot(xt, ltmap, '--', color=cols[0], linewidth=2) )
        if truemodel is not None:
            h.append( plt.plot(truemodel.x, truemodel.l, 'k-') )

        if (mapgp is not None) and (truemodel is not None):
            plt.legend(h, ['Samples','MAP lengthscale', 'True lengthscale'], loc='upper left')
        elif truemodel is not None:
            plt.legend(h, ['Samples','True lengthscale'], loc='upper left')
        elif mapgp is not None:
            plt.legend(h, ['Samples','MAP lengthscale'], loc='upper left')
        else:
            plt.legend(h, ['Samples'])

        plt.ylabel('value')
        plt.ylim(0, np.quantile(np.max(lts, axis=1),0.99))
        plt.title('Lengthscale posterior')

        plt.subplot(squares,1,4)
        h = [None]
        for i in range(n):
            h[0] = plt.fill(np.append(xt, np.flip(xt)), np.append(ots[i] + lw_o, np.flip(ots[i] - lw_o)), color=cols[2], alpha=alpha_o/n, edgecolor=None)
        if mapgp is not None:
            h.append( plt.plot(xt, otmap, '--', color=cols[2], linewidth=2) )
        if truemodel is not None:
            h.append( plt.plot(truemodel.x, truemodel.o, 'k-') )

        if (mapgp is not None) and (truemodel is not None):
            plt.legend(h, ['Samples','MAP noise variance', 'True noise variance'], loc='upper left')
        elif truemodel is not None:
            plt.legend(h, ['Samples','True noise variance'], loc='upper left')
        elif mapgp is not None:
            plt.legend(h, ['Samples','MAP noise variance'], loc='upper left')
        else:
            plt.legend(h, ['Samples'])

#        plt.xlabel('time')
        plt.ylabel('value')
        plt.ylim(0, np.quantile(np.max(ots, axis=1),0.99))
        plt.title('Noise variance posterior')

        plt.subplot(squares,1,5)
        h = [None]
        for i in range(n):
            h[0] = plt.fill(np.append(xt, np.flip(xt)), np.append(sts[i] + lw_s, np.flip(sts[i] - lw_s)), color=cols[1], alpha=alpha_s/n, edgecolor=None)
        if mapgp is not None:
            h.append( plt.plot(xt, stmap, '--', color=cols[1], linewidth=2) )
        if truemodel is not None:
            h.append( plt.plot(truemodel.x, truemodel.s, 'k-') )

        if (mapgp is not None) and (truemodel is not None):
            plt.legend(h, ['Samples','MAP signal variance', 'True variance'], loc='upper left')
        elif truemodel is not None:
            plt.legend(h, ['Samples','True signal variance'], loc='upper left')
        elif mapgp is not None:
            plt.legend(h, ['Samples','MAP signal variance'], loc='upper left')
        else:
            plt.legend(h, ['Samples'])

        plt.ylabel('value')
        plt.xlabel('time')
        plt.ylim(0, np.quantile(np.max(sts, axis=1),0.99))
        plt.title('Signal variance posterior')


#//#
#//# testerrors.m
#//#

def testerrors(gp, samples=np.empty(0)):

    m = samples.shape[0]
#    yts = pars.yts * pars.yscale + pars.ybias
#    ytr = pars.ytr * pars.yscale + pars.ybias

    nmse = 0

    if np.all(samples):
        mse = np.zeros((m,1))
#        nmse = np.zeros((m,1))
        nlpd = np.zeros((m,1))

        for i in range(m):
            gp = parsesample(samples[i], gp)
            mts, stds, _, _, ots = nsgpposterior(gp, gp.xts)

            # compute over samples
            mse[i] = computemse(mts, yts)
#            nmse[i] = computenmse(mts, yts, ytr)
            nlpd[i] = computenlpd(mts, yts, stds**2 + ots**2)
    else:
        # compute posterior for test points
        ft, stds, _, _, ot = nsgpposterior(gp, gp.xts)

        # compute statistics: MLL, NMSE, NLPD
        mse = computemse(ft, gp.yts)
#        nmse = computenmse(ft, gp.yts, gp.ytr)
        nlpd = computenlpd(ft, gp.yts, stds**2 + ot**2)

    return mse, nmse, nlpd
