from data_utils import read_graphfile
import numpy as np
import pandas as pd
import os
import sys

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold
from sklearn.svm import SVC, LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
from sklearn.metrics import accuracy_score
from sklearn.manifold import TSNE

def evaluate_embedding(embeddings, labels):

    labels = preprocessing.LabelEncoder().fit_transform(labels)
    x, y = np.array(embeddings), np.array(labels)
    print(x.shape, y.shape)
    
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=None)
    accuracies = []
    for train_index, test_index in kf.split(x, y):

        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        search=True
        if search:
            params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
            classifier = GridSearchCV(SVC(), params, cv=5, scoring='accuracy', verbose=0)
        else:
            classifier = SVC(C=10)
        classifier.fit(x_train, y_train)
        accuracies.append(accuracy_score(y_test, classifier.predict(x_test)))

    svm_accuracies = np.mean(accuracies)

    accuracies = []
    for train_index, test_index in kf.split(x, y):

        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        search=True
        if search:
            params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
            classifier = GridSearchCV(LinearSVC(), params, cv=5, scoring='accuracy', verbose=0)
        else:
            classifier = LinearSVC(C=10)
        classifier.fit(x_train, y_train)
        accuracies.append(accuracy_score(y_test, classifier.predict(x_test)))
    print('LinearSvc', np.mean(accuracies))
    print('svc', svm_accuracies)

# def get_mutag():
    # emb = []
    # with open('data/results/output.txt', 'r') as f:
        # for line in f:
            # emb.append(list(map(float, [x for x in line.strip().split()])))

    # ret = []
    # for i in range(188):
        # with open('./data/mutag/mutag_{}.graph'.format(i+1), 'r') as f:
            # x = f.readlines()
        # ret.append(int(x[-1].strip()))
    # return emb, ret


if __name__ == '__main__':
    # x, y = get_mutag()
    emb = []
    with open('data/results/{}_output.txt'.format(sys.argv[1]), 'r') as f:
        for line in f:
            emb.append(list(map(float, [x for x in line.strip().split()])))

    with open('../data/{}_label.txt'.format(sys.argv[1]), 'r') as f:
        y = f.readlines()
    y = [int(x.strip()) for x in y]
    
    evaluate_embedding(emb, y)
