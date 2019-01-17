import numpy as np
import networkx as nx
from glob import glob
from tqdm import tqdm
import os
import subprocess
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.svm import SVC, LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
from sklearn.metrics import accuracy_score

def load_data(ds_name, use_node_labels):
    node2graph = {}
    Gs = []

    with open("../data/%s/%s_graph_indicator.txt"%(ds_name,ds_name), "r") as f:
            c = 1
            for line in f:
                    node2graph[c] = int(line[:-1])
                    if not node2graph[c] == len(Gs):
                            Gs.append(nx.Graph())
                    Gs[-1].add_node(c)
                    c += 1

    with open("../data/%s/%s_A.txt"%(ds_name,ds_name), "r") as f:
            for line in f:
                    edge = line[:-1].split(",")
                    edge[1] = edge[1].replace(" ", "")
                    Gs[node2graph[int(edge[0])]-1].add_edge(int(edge[0]), int(edge[1]))

    if use_node_labels:
            with open("../data/%s/%s_node_labels.txt"%(ds_name,ds_name), "r") as f:
                    c = 1
                    for line in f:
                            node_label = int(line[:-1])
                            Gs[node2graph[c]-1].node[c]['label'] = node_label
                            c += 1

    labels = []
    with open("../data/%s/%s_graph_labels.txt"%(ds_name,ds_name), "r") as f:
            for line in f:
                    labels.append(int(line[:-1]))

    labels  = np.array(labels, dtype = np.float)
    return Gs, labels

def preprocess(DS):
    Gs, labels = load_data(DS, False)
    print('number of graphs', len(Gs))

    datadir = '../data/{}'.format(DS)
    try:
        os.mkdir(datadir)
    except Exception as e:
        print(e)

    assert len(Gs) == len(labels)
    f = open('../data/{}.Labels'.format(DS), 'w')
    for graphidx, G in tqdm(enumerate(Gs)):
        nx.write_gexf(G, '{}/{}.gexf'.format(datadir, graphidx))
        f.write('{}.gexf {}\n'.format(graphidx, int(labels[graphidx])))
    f.close()



if __name__ == '__main__':
    import sys
    preprocess(sys.argv[1])
    # preprocess('ENZYMES')
    # preprocess('DD')
    # preprocess('REDDIT-BINARY')
    # preprocess('COLLAB')
    # preprocess('REDDIT-MULTI-5K')
    # preprocess('IMDB-BINARY')
    # preprocess('IMDB-MULTI')
