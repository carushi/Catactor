import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/../LassoVariants/AlternateLasso')

# sys.path.append('../../LassoVariants/AlternateLasso')
import pickle
import numpy as np
from AlternateLinearModel import AlternateLogisticLasso
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
import pandas as pd
from copy import deepcopy
from sklearn.metrics import roc_auc_score
from multiprocessing import Pool


SEED = 42

def apply_alternate_LL(X, y, feature_names=[], rho=0.1, tol=1e-4, save='test.npy'):
    mdl = AlternateLogisticLasso(rho=rho, tol=tol, verbose=True, save=save)
    mdl = mdl.fit(X, y, feature_names)
    #print(mdl.a_, mdl.lasso_.coef_)
    return mdl

def parse_mdl_predictor(mdl):
    if len(mdl.alternates_) == 0:
        for i, coef in enumerate(mdl.a_):
            yield i, i, 0, 0, coef
    else:
        for i, best_features in enumerate(mdl.alternates_):
            fname, coef = best_features[0], best_features[1]
            yield  fname, i, 0, 0, coef
            for j, semi in enumerate(best_features[2]):
                sname, sscore, scoef = semi[0], semi[1], semi[2]
                yield sname, i, j+1, sscore-mdl.obj_, scoef
    yield 'intercept', -1, 0, 0, mdl.b_

METHODS = ['la', 'logistic', 'rf', 'svm']
# METHODS = ['la', 'logistic', 'rf']

def get_classifier(classifier, kwargs={}):
    mdl = None
    if classifier == 'la' or classifier == 'laobest':
        mdl = AlternateLogisticLasso(**kwargs)
    if classifier == 'logistic':
        mdl = LogisticRegression(penalty='l1', solver='liblinear')
    elif classifier == 'rf':
        mdl = RandomForestClassifier(n_estimators=100)
    elif classifier == 'svm':
        mdl = SVC(kernel='rbf', degree=3, probability=True, random_state=SEED, gamma='auto')
    else:
        pass
    return mdl

def make_other_clfs_each(X, c, clusters, classifier, output):
    y = np.array([1 if cl == c else 0 for cl in clusters])
    mdl = get_classifier(classifier)
    mdl = mdl.fit(X, y)
    print('Best auc', output, c, classifier, roc_auc_score(y, mdl.predict_proba(X)[:,1]))
    return (c, deepcopy(mdl))

def make_other_clfs_parallel(X, clusters, output, classifier='logistic', cores=1): # logistic, svm, random forest
    global SEED
    pool = Pool(processes=cores)
    print(cores, classifier)
    # X = X[0:100,:]
    # clusters = clusters[0:100]
    cluster_list = sorted(list(set(clusters)))
    for i in range(0, len(cluster_list), cores):
        print(str(i)+'-'+str(min(len(cluster_list), i+cores))+' clusters')
        sys.stdout.flush()
        result = pool.starmap_async(make_other_clfs_each, [(deepcopy(X), c, clusters, classifier, output) for c in cluster_list[i:min(len(cluster_list), i+cores)]])
        clfs = {}
        for (c, mdl) in result.get():
            clfs[c] = mdl
    with open(output+'_'+classifier+'.npy', 'wb') as f:
        pickle.dump(clfs, f)
    pool.close()
    pool.join()
    del pool

def make_other_clfs(X, clusters, output, classifier='logistic', cores=1, saved=False): # logistic, svm, random forest
    global SEED
    if cores > 1:
        make_other_clfs_parallel(X, clusters, output, classifier, cores)
        return
    clfs = {}
    ofname = output+'_'+classifier+'.npy'
    if (not saved) or (not os.path.exists(ofname)):
        for c in set(clusters):
            y = np.array([1 if cl == c else 0 for cl in clusters])
            mdl = get_classifier(classifier)
            mdl = mdl.fit(X, y)
            print('Best auc', output, c, classifier, roc_auc_score(y, mdl.predict_proba(X)[:,1]))
            clfs[c] = deepcopy(mdl)
        with open(ofname, 'wb') as f:
            pickle.dump(clfs, f)

def read_other_clfs(output, classifier='logistic'):
    if os.path.exists(output+'_'+classifier+'.npy'):
        with open(output+'_'+classifier+'.npy', 'rb') as f:
            return pickle.load(f)
    return None

def read_lasso_clf(header):
    if os.path.exists(header+'.npy'):
        with open(header+'.npy', 'rb') as f:
            return pickle.load(f)
    return None

def make_lasso_matrix_each(X, c, clusters, feature_names):
    print('make_lasso_each', c, X.shape)
    kwargs = {'rho':1e-5, 'tol':1e-4, 'check':0.90, 'max_alternates':-1, 'maxitr':1000}
    mdl = get_classifier('la', kwargs)
    y = np.array([True if cl == c else False for cl in clusters])
    mdl = mdl.fit(X, y, featurename=feature_names)
    # mdl.predict(X)
    temp = pd.DataFrame([list(x) for x in parse_mdl_predictor(mdl)], columns=['index', 'dimension', 'order', 'score', 'coef'])
    temp = temp.assign(gene=[feature_names[int(x)] if x != 'intercept' else x for x in temp.loc[:,'index']])
    temp = temp.assign(celltype=c)
    print(c, mdl, temp)
    return (c, deepcopy(mdl), temp)

def make_lasso_matrix_parallel(X, clusters, header, feature_names=[], cores=1):
    print(cores, 'la')
    print(X.shape, clusters[0:10], len(set(clusters)))
    pmatrix = None
    clfs = {}
    pool = Pool(processes=cores)
    cluster_list = sorted(list(set(clusters)))
    for i in range(0, len(cluster_list), cores):
        print(str(i)+'-'+str(min(len(cluster_list), i+cores))+' clusters')
        sys.stdout.flush()
        result = pool.starmap_async(make_lasso_matrix_each, [(deepcopy(X), c, clusters, feature_names) for c in cluster_list[i:min(len(cluster_list), i+cores)]])
        for (c, mdl, temp) in result.get():
            if pmatrix is None:
                pmatrix = temp
            else:
                pmatrix = pd.concat([pmatrix, temp], axis=0, ignore_index=True)
            clfs[c] = mdl
    with open(header+'.npy', 'wb') as f:
        pickle.dump(clfs, f)
    if pmatrix is not None:
        pmatrix.to_csv(header+'.csv')
    pool.close()
    pool.join()
    del pool

def make_lasso_matrix(X, clusters, header, feature_names=[], cores=1, saved=False):
    if cores > 1:
        make_lasso_matrix_parallel(X, clusters, header, feature_names, cores)
        return
    pmatrix = None
    clfs = {}
    kwargs = {'rho':1e-5, 'tol':1e-4, 'check':0.90, 'max_alternates':-1, 'maxitr':1000}
    ofname = header+'.npy'
    if (not saved) or (not os.path.exists(ofname)):
        for c in sorted(list(set(clusters))):
            mdl = get_classifier('la', kwargs)
            y = np.array([True if cl == c else False for cl in clusters])
            mdl = mdl.fit(X, y, featurename=feature_names)
            mdl.predict(X)
            # mdl = apply_alternate_LL(X, y, feature_names, save='')
            temp = pd.DataFrame([list(x) for x in parse_mdl_predictor(mdl)], columns=['index', 'dimension', 'order', 'score', 'coef'])
            temp = temp.assign(gene=[feature_names[int(x)] if x != 'intercept' else x for x in temp.loc[:,'index']])
            temp = temp.assign(celltype=c)
            if pmatrix is None:
                pmatrix = temp
            else:
                pmatrix = pd.concat([pmatrix, temp], axis=0, ignore_index=True)
            clfs[c] = deepcopy(mdl)
        with open(ofname, 'wb') as f:
            pickle.dump(clfs, f)
        if pmatrix is not None:
            pmatrix.to_csv(header+'.csv')

if __name__ == "__main__":
    seed = 0
    num = 1000
    dim = 2
    dim_extra = 2

    np.random.seed(seed)
    X = np.random.randn(num, dim + dim_extra)
    for i in range(dim_extra):
        X[:, dim + i] = X[:, 0] + 0.5 * np.random.randn(num)
    y = X[:, 0] + 0.3 * X[:, 1] + 0.5 * np.random.randn(num) > 0
    print(apply_alternate_LL(X, y).predict(X, 1))
    make_lasso_matrix(X, y, 'temp_mat.csv', [0, 1, 2, 3])


