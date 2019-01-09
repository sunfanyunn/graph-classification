from glob import glob
from grakel import Graph
from grakel import GraphKernel
from grakel import datasets
from graph_tool import load_graph
from matplotlib import pylab as pl
from sklearn import svm
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.model_selection import train_test_split
from time import time
from tqdm import tqdm
import numpy as np
import sys


def get_data(graphs):

    data = []
    for idx, G in tqdm(enumerate(graphs)):
        
        D = []
        D.append({ (e[0]+1, e[1]+1) for e in G.edges()})

        node_labels = {}
        for u in G.nodes():
            try:
                n = op_map[G.node[u]['label']]
            except:
                n = op_map[G.node[u]['label']] = len(op_map)

            node_labels[int(u)+1] = n

        D.append(node_labels)
        
        edge_labels = { (e[0]+1, e[1]+1): 1 for e in G.edges()}
        D.append(edge_labels)

        D = Graph(D[0], node_labels=node_labels, edge_labels=None, graph_format='auto')
        data.append(D)


    return data

if __name__ == '__main__':
    op_map = {}

    import sys
    from utils import read_graphfile
    G, y = read_graphfile('../data/', sys.argv[1])
    method = sys.argv[2]
    unique, counts = np.unique(y, return_counts=True)
    print(unique, counts)
    # input()
    G = get_data(G)
    
    if method == 'wl':
        for niter in [1,2,3,4,5]:
            ultimate_accs = []
            for _ in range(5):
                kf = StratifiedKFold(n_splits=10, shuffle=True)
                accs = []
                for train_index, test_index in kf.split(G, y):
            
                    start = time()
            
                    G_train = [G[idx] for idx in train_index]
                    y_train = [y[idx] for idx in train_index]
                    G_test = [G[idx] for idx in test_index]
                    y_test = [y[idx] for idx in test_index]
            
                    # Initialise a weifeiler kernel, with a dirac base_kernel.
                    gk = GraphKernel(kernel=[{"name": "weisfeiler_lehman", "niter": niter},
                                             {"name": "subtree_wl"}], normalize=True)
                    
                    # Calculate the kernel matrix.
                    K_train = gk.fit_transform(G_train)
                    K_test = gk.transform(G_test)
                    
                    # Initialise an SVM and fit.
                    clf = svm.SVC(kernel='precomputed', C=1)
                    params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
                    clf = GridSearchCV(svm.SVC(kernel='precomputed'), params, cv=10, scoring='accuracy', verbose=0)
                    # print('fitting ...')
                    clf.fit(K_train, y_train)
                    
                    # Predict and test.
                    # print('predicting ...')
                    y_pred = clf.predict(K_test)
                    # y_pred = np.random.randint(0, 3, len(y_test))
                    
                    # Calculate accuracy of classification.
                    acc = accuracy_score(y_test, y_pred)
                    accs.append(acc)
                    
                    end = time()
                    # print("Accuracy:", str(round(acc*100, 2)), "% | Took:",
                          # str(round(end - start, 2)), "s")
            
                print('wl', niter, np.mean(accs), np.std(accs))
                ultimate_accs.append(np.mean(accs))
            print(np.mean(ultimate_accs), np.std(ultimate_accs))
    
    if method == 'graphlet':
        name = "graphlet_sampling"
        ultimate_accs = []
        for n_sample in [1, 2, 3, 4, 5]:
            ultimate_accs = []
            for _ in range(5):
                accs = []
                kf = StratifiedKFold(n_splits=10, shuffle=True)
                for train_index, test_index in kf.split(G, y):
            
                    start = time()
            
                    G_train = [G[idx] for idx in train_index]
                    y_train = [y[idx] for idx in train_index]
                    G_test = [G[idx] for idx in test_index]
                    y_test = [y[idx] for idx in test_index]
            
                    gk = GraphKernel(kernel=[{"name": name, "sampling":{'n_samples':n_sample}}], normalize=True)
                    
                    # Calculate the kernel matrix.
                    K_train = gk.fit_transform(G_train)
                    K_test = gk.transform(G_test)
                    
                    # Initialise an SVM and fit.
                    # clf = svm.SVC(kernel='precomputed', C=1)
                    params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
                    clf = GridSearchCV(svm.SVC(kernel='precomputed'), params, cv=10, scoring='accuracy', verbose=0)
                    clf.fit(K_train, y_train)
                    
                    # Predict and test.
                    y_pred = clf.predict(K_test)
                    
                    # Calculate accuracy of classification.
                    acc = accuracy_score(y_test, y_pred)
                    accs.append(acc)
                    
                    end = time()
                    # print("Accuracy:", str(round(acc*100, 2)), "% | Took:",
                          # str(round(end - start, 2)), "s")
            
                print(name, np.mean(accs), np.std(accs))
                ultimate_accs.append(np.mean(accs))
            print(np.mean(ultimate_accs), np.std(ultimate_accs))

    if method == 'shortest':
        ultimate_accs = []
        for _ in range(5):
            name = "shortest_path"
            kf = StratifiedKFold(n_splits=10, shuffle=True)
            accs = []
            for train_index, test_index in kf.split(G, y):
        
                start = time()
        
                G_train = [G[idx] for idx in train_index]
                y_train = [y[idx] for idx in train_index]
                G_test = [G[idx] for idx in test_index]
                y_test = [y[idx] for idx in test_index]
        
                gk = GraphKernel(kernel=[{"name": name}], normalize=True)
                
                # Calculate the kernel matrix.
                K_train = gk.fit_transform(G_train)
                K_test = gk.transform(G_test)
                
                # Initialise an SVM and fit.
                # clf = svm.SVC(kernel='precomputed', C=1)
                params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
                clf = GridSearchCV(svm.SVC(kernel='precomputed'), params, cv=10, scoring='accuracy', verbose=0)
                clf.fit(K_train, y_train)
                
                # Predict and test.
                y_pred = clf.predict(K_test)
                
                # Calculate accuracy of classification.
                acc = accuracy_score(y_test, y_pred)
                accs.append(acc)
                
                end = time()
                # print("Accuracy:", str(round(acc*100, 2)), "% | Took:",
                      # str(round(end - start, 2)), "s")
        
            # print(name, np.mean(accs), np.std(accs))
            ultimate_accs.append(np.mean(accs))
        print(np.mean(ultimate_accs), np.std(ultimate_accs))

    if method == 'walk':
        ultimate_accs = []
        for _ in range(5):
            print("random_walk")
            kf = StratifiedKFold(n_splits=10, shuffle=True)
            accs = []
            for train_index, test_index in kf.split(G, y):
        
                start = time()
        
                G_train = [G[idx] for idx in train_index]
                y_train = [y[idx] for idx in train_index]
                G_test = [G[idx] for idx in test_index]
                y_test = [y[idx] for idx in test_index]
        
                gk = GraphKernel(kernel=[{"name": "random_walk"}], normalize=True)
                
                # Calculate the kernel matrix.
                K_train = gk.fit_transform(G_train)
                K_test = gk.transform(G_test)
                
                # Initialise an SVM and fit.
                # clf = svm.SVC(kernel='precomputed', C=1)
                params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
                clf = GridSearchCV(svm.SVC(kernel='precomputed'), params, cv=10, scoring='accuracy', verbose=0)
                clf.fit(K_train, y_train)
                
                # Predict and test.
                y_pred = clf.predict(K_test)
                
                # Calculate accuracy of classification.
                acc = accuracy_score(y_test, y_pred)
                accs.append(acc)
                
                end = time()
                # print("Accuracy:", str(round(acc*100, 2)), "% | Took:",
                      # str(round(end - start, 2)), "s")
        
            print(name, np.mean(accs), np.std(accs))
            ultimate_accs.append(np.mean(accs))
        print(np.mean(ultimate_accs), np.std(ultimate_accs))
