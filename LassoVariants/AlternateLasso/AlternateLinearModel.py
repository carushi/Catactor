# -*- coding: utf-8 -*-
"""
@author: Satoshi Hara

(Class)
> AlternateLasso(maxitr=1000, tol=1e-4, rho=0.1, verbose=False, save='') # for linear regression
> AlternateLogisticLasso(maxitr=1000, tol=1e-4, rho=0.1, verbose=False, save='') # for linear binary classification
    maxitr      : maximum number of iterations for Lasso
    tol         : tolerance parameter to stop the iterative optimization
    rho         : regularization parameter
    verbose     : print the feature search process when 'verbose=True'
    save        : save the model during the feature search into the specified file, e.g., save='model.npy'

(Method)
> AlternateLasso.fit(X, y, featurename=[])
> AlternateLogisticLasso.fit(X, y, featurename=[])
    X           : numpy array of size num x dim
    y           : numpy array of size num
    featurename : name of features
"""

import numpy as np
from sklearn.linear_model import Lasso
from sklearn.linear_model import LogisticRegression
from sklearn.externals import joblib

class AlternateLasso(object):
    def __init__(self, maxitr=1000, tol=1e-4, rho=0.1, verbose=False, save=''):
        self.maxitr_ = maxitr
        self.tol_ = tol
        self.rho_ = rho
        self.verbose_ = verbose
        self.save_ = save
        self.lasso_ = Lasso(max_iter=self.maxitr_, tol=self.tol_, alpha=self.rho_)
        
    def __str__(self):
        s = ''
        for a in self.alternates_:
            s += 'Feature: %s, Coef. = %f\n' % (self.featurename_[a[0]], self.a_[a[0]])
            if len(a[2]) == 0:
                s += '\t ** No Alternate Features **\n'
                continue
            for b in a[2]:
                s += '\t Alternate Feature: %s, Score = %f, Coef. = %f\n' % (self.featurename_[b[0]], b[2] - self.obj_, b[1])
        return s
        
    def fit(self, x, y, featurename=[]):
        self.dim_ = x.shape[1]
        self.setfeaturename(featurename)
        y = np.array(y).reshape((1, y.shape[0]))        
        # lasso
        self.lasso_.fit(x, y)
        self.a_ = self.lasso_.coef_
        self.b_ = self.lasso_.intercept_
        self.obj_ = self.__getObj(x, y, self.a_, self.b_)
        nonzeros = self.a_.nonzero()[0]
        
        # alternate lasso
        self.alternates_ = []
        p = self.__getP(x, y, self.a_, self.b_)
        q = self.__getQ(x)
        if self.verbose_:
            print('> [feature name, # of alternate feature candidates]')
        for i, d in enumerate(nonzeros):
            aa = self.a_.copy()
            aa[d] = 0
            obj0 = self.__getObj(x, y, aa, self.b_)
            pp = p.copy()
            pp -= self.__getG(x, d) * self.a_[d]
            r = pp / q
            da = np.sign(- r) * np.maximum(0, np.abs(r) - self.rho_ / q)
            da[nonzeros] = 0
            da[np.isnan(da)] = 0
            gg = []
            if self.verbose_:
                print('> [%s, %d]' % (self.featurename_[d], da.nonzero()[0].size))
            for dd in da.nonzero()[0]:
                aa[dd] = da[dd]
                objt = self.__getObj(x, y, aa, self.b_)
                aa[dd] = 0
                gg.append((dd, da[dd], objt))
            self.alternates_.append((d, obj0, gg))
            if len(self.save_) > 0:
                joblib.dump(self, self.save_, compress=9)
    
    def predict(self, x, idx=-1):
        if idx == -1:
            return x.dot(self.a_) + self.b_
        else:
            a = self.a_.copy()
            b = self.b_
            g = self.g_[idx]
            a[g[0]] = 0
            z = np.zeros((x.shape[0], len(g[2])+1))
            z[:, 0] = x.dot(a) + b >= 0
            for i, gg in enumerate(g[2]):
                a[gg[0]] = gg[1]
                z[:, i+1] = x.dot(a) + b
                a[gg[0]] = 0
            return z
    
    def evaluate(self, x, y):
        z = self.predict(x)
        mse = [np.mean((y - z)**2)]
        for i in range(len(self.alternates_)):
            z = self.predict(x, idx=i)
            mse.append(list(np.mean((y[:, np.newaxis] - z)**2, axis=0)))
        return mse
    
    def setfeaturename(self, featurename):
        if len(featurename) > 0:
            self.featurename_ = featurename
        else:
            self.featurename_ = []
            for d in range(self.dim_):
                self.featurename_.append('x_%d' % (d+1,))
    
    def __getObj(self, x, y, a, b):
        z = x.dot(a) + b
        f = np.mean((y - z)**2) / 2
        return f + self.rho_ * np.sum(np.abs(a))
    
    def __getG(self, x, idx):
        num = x.shape[0]
        return x.T.dot(x[:, idx]) / num
    
    def __getH(self, x, y, b):
        num = y.size
        return - x.T.dot(y - b) / num
    
    def __getP(self, x, y, a, b):
        p = self.__getH(x, y, b)
        for d in a.nonzero()[0]:
            g = self.__getG(x, d)
            p += g * a[d]
        return p
    
    def __getQ(self, x):
        return np.mean(x**2, axis=0)
        
class AlternateLogisticLasso(object):
    def __init__(self, maxitr=1000, tol=1e-4, rho=0.1, lasso_C=1., lasso_tol=1e-4, lasso_itr=100, check=0.99, max_alternates=-1, verbose=False, save=''):
        self.maxitr_ = maxitr
        self.tol_ = tol
        self.rho_ = rho
        self.check_ = check
        self.verbose_ = verbose
        self.save_ = save
        #self.lasso_ = LogisticRegression(max_iter=self.maxitr_, tol=self.tol_, C=1/self.rho_, penalty='l1', solver='liblinear')
        self.lasso_ = LogisticRegression(penalty='l1', solver='liblinear')
        self.max_alternates = max_alternates
        
    def __str__(self):
        s = ''
        for a in self.alternates_:
            s += 'Feature: %s, Coef. = %f\n' % (self.featurename_[a[0]], self.a_[a[0]])
            if len(a[2]) == 0:
                s += '\t ** No Alternate Features **\n'
                continue
            for b in a[2]:
                s += '\t Alternate Feature: %s, Score = %f, Coef. = %f\n' % (self.featurename_[b[0]], b[2] - self.obj_, b[1])
        return s
        
    def fit(self, x, y, featurename=[]):
        self.dim_ = x.shape[1]
        self.setfeaturename(featurename)
        # lasso
        self.lasso_.C = 1 / (self.rho_ * x.shape[0])
        self.lasso_.fit(x, y)
        self.a_ = self.lasso_.coef_[0, :]
        self.b_ = self.lasso_.intercept_[0]
        self.obj_ = self.__getObj(x, y, self.a_, self.b_)
        nonzeros = self.a_.nonzero()[0]
        
        # alternate lasso
        self.alternates_ = []
        self.count_ = 0
        if self.verbose_:
            print('> [feature name, # of alternate feature candidates]')
        for i, d in enumerate(nonzeros):
            aa = self.a_.copy()
            aa[d] = 0
            obj0 = self.__getObj(x, y, aa, self.b_)
            p = self.__getG(x, y, aa, self.b_)
            p[nonzeros] = 0
            p = (np.abs(p) / self.rho_ > self.check_)
            # print(np.abs(p)/self.rho_, self.check_)
            gg = []
            if self.verbose_:
                print('> [%s, %d]' % (self.featurename_[d], p.nonzero()[0].size))
            for dd in p.nonzero()[0]:
                self.count_ += 1
                da = self.__fitLogReg(x, y, aa, self.b_, dd, maxitr=5000)
                if da == 0:
                    continue
                aa[dd] = da
                objt = self.__getObj(x, y, aa, self.b_)
                aa[dd] = 0
                gg.append((dd, da, objt))
            self.alternates_.append((d, obj0, gg))
            if len(self.save_) > 0:
                joblib.dump(self, self.save_, compress=9)
        return self
        
    def predict(self, x, idx=-1, max_alternates=-1):
        if idx == -1:
            return x.dot(self.a_) + self.b_ >= 0
        return self.predict_proba_mean(x, idx, max_alternates) >= 0

    def predict_proba_mean(self, x, idx=-1, max_alternates=-1):
        if idx == -1:
            max_idx = len(self.alternates_)
        else:
            max_idx = min(len(self.alternates_), idx)
        result = []
        #print('alternates', max_idx, self.a_.max(axis=-1))
        Z = np.array([self.predict_proba(x, idx=i, max_alternates=max_alternates) for i in range(-1, max_idx)])
        #print('Z', Z.shape, Z)
        return Z.mean(axis=0)

    def predict_proba(self, x, idx=-1, max_alternates=-1):
        if idx == -1:
            return x.dot(self.a_) + self.b_
        else:
            g = self.alternates_[idx]
            if x[:,g[0]].sum() > 0 or len(g[2]) == 0:
                return x.dot(self.a_)+self.b_
            print('Not detected', self.featurename_[idx], 'dimension', idx)
            a = self.a_.copy()
            b = self.b_
            a[g[0]] = 0
            limits = [(self.max_alternates if max_alternates < 0 else max_alternates), len(g[2])]
            max_it = min([x for x in limits if x >= 0])
            z = np.zeros((x.shape[0], len(g[2])+1))
            z[:, 0] = x.dot(a) + b
            score_list = [(i, gg[2]) for i, gg in enumerate(g[2])]
            score_ordered = sorted(score_list, key=lambda x: x[1], reverse=True)
            observed = []
            for n, (i, score) in enumerate(score_ordered):
                if n >= max_it: break
                gg = g[2][i]
                if x[:,gg[0]].sum() == 0 and idx >= 0:
                    continue
                a[gg[0]] = gg[1]
                z[:, n+1] = x.dot(a) + b
                a[gg[0]] = 0
                observed.append(n)
            print(len(observed), observed)
            if len(observed) == 0:
                return x.dot(self.a_)+self.b_
            return z[:, np.array(observed)].mean(axis=1)

    def predict_proba_mean_weight(self, x, idx=-1, max_alternates=5):
        if idx == -1:
            return x.dot(self.a_) + self.b_
        else:
            g = self.alternates_[idx]
            if x[:,g[0]].sum() > 0 or len(g[2]) == 0:
                return x.dot(self.a_)+self.b_
            print('Not detected', self.featurename_[idx], 'dimension', idx)
            a = self.a_.copy()
            b = self.b_
            a[g[0]] = 0
            limits = [(self.max_alternates if max_alternates < 0 else max_alternates), len(g[2])]
            max_it = min([x for x in limits if x >= 0])
            z = np.zeros((x.shape[0], len(g[2])+1))
            z[:, 0] = (x.dot(a) + b)/2.
            score_list = [(i, gg[2]) for i, gg in enumerate(g[2])]
            score_ordered = sorted(score_list, key=lambda x: x[1], reverse=True)
            observed = []
            for n, (i, score) in enumerate(score_ordered):
                if n >= max_it: break
                gg = g[2][i]
                if x[:,gg[0]].sum() == 0 and idx >= 0:
                    continue
                a[gg[0]] = gg[1]
                z[:, n+1] = (x.dot(a) + b)/max_it
                a[gg[0]] = 0
                observed.append(n)
            print(len(observed), observed)
            z[:, np.array(observed)] /= len(observed)*2
            observed.append(0)
            if len(observed) == 0:
                return x.dot(self.a_)+self.b_
            return z[:, np.array(observed)].sum(axis=1)

    def evaluate(self, x, y):
        z = self.predict(x)
        mcr = [1 - np.mean(y == z)]
        for i in range(len(self.g_)):
            z = self.predict(x, idx=i)
            mcr.append(list(1 - np.mean(y[:, np.newaxis] == z, axis=0)))
        return mcr
    
    def setfeaturename(self, featurename):
        if len(featurename) > 0:
            self.featurename_ = featurename
        else:
            self.featurename_ = []
            for d in range(self.dim_):
                self.featurename_.append('x_%d' % (d+1,))
    
    def __getObj(self, x, y, a, b):
        z = x.dot(a) + b
        f = np.mean(np.log(1 + np.exp(- (2 * y - 1) * z)))
        return f + self.rho_ * np.sum(np.abs(a))
    
    def __getG(self, x, y, a, b):
        num = y.size
        z = x.dot(a) + b
        r = np.exp(- (2*y - 1) * z)
        return - x.T.dot((2 * y - 1) * r / (1 + r)) / num
        
    def __fitLogReg(self, x, y, a, b, idx, maxitr=1000, tol=1e-10):
        z = np.squeeze(np.array(x.dot(a) + b))
        alpha = (2 * y - 1) * x[:, idx]
        beta = (2 * y - 1) * z
        c = 0
        obj = self.__objLogReg(alpha, beta, c)
        smax = 10
        r = 0.3
        rc = 10
        for itr in range(self.maxitr_):
            dc = self.__getLogGrad(alpha, beta, c)
            s = smax
            for count in range(rc):
                cc = c - s * dc
                objt = self.__objLogReg(alpha, beta, cc)
                if objt < obj:
                    break
                s *= r
            cnew = c - s * dc
            cnew = np.sign(cnew) * max(0, np.abs(cnew) - s * self.rho_)
            objnew = self.__objLogReg(alpha, beta, cnew)
            if obj - objnew < tol:
                if objnew < obj:
                    c = cnew
                break
            c = cnew
            obj = objnew
        return c
        
    def __objLogReg(self, alpha, beta, c):
        f = np.mean(np.log(1 + np.exp(- beta - alpha * c)))
        return f + self.rho_ * np.abs(c)
    
    def __getLogGrad(self, alpha, beta, c):
        num = alpha.size
        r = np.exp(- beta - alpha * c)
        return - alpha.T.dot(r / (1 + r)) / num
        
