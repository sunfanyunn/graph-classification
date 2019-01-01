"""
Please visit www.mit.edu/~pinary/kdd for more information.

The scripts is licensed with MIT License, which you can find a copy below:

Copyright (c) 2016 Pinar Yanardag

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
import numpy as np
import multiprocessing as mp
import networkx as nx
from gensim.models import Word2Vec
from itertools import chain, combinations
from collections import defaultdict
import sys, copy, time, math, pickle
import itertools
import scipy.io
import pynauty
import random
from scipy.spatial.distance import pdist, squareform

#random.seed(314124)
#np.random.seed(2312312)

import networkx as nx
import numpy as np
from tqdm import tqdm
from glob import glob
import re
import os

def load_data(dataname, datadir='../data', max_nodes=None):
    ''' Read data from https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
        graph index starts with 1 in file
    Returns:
        List of networkx objects with graph and node labels
    '''
    prefix = os.path.join(datadir, dataname, dataname)
    filename_graph_indic = prefix + '_graph_indicator.txt'
    # index of graphs that a given node belongs to
    graph_indic={}
    with open(filename_graph_indic) as f:
        i=1
        for line in f:
            line=line.strip("\n")
            graph_indic[i]=int(line)
            i+=1

    filename_nodes=prefix + '_node_labels.txt'
    node_labels=[]
    try:
        with open(filename_nodes) as f:
            for line in f:
                line=line.strip("\n")
                node_labels+=[int(line)]
        # node_labels = LabelEncoder().fit_transform(node_labels)
    except IOError:
        print('No node labels')

    filename_node_attrs=prefix + '_node_attributes.txt'
    node_attrs=[]
    try:
        with open(filename_node_attrs) as f:
            for line in f:
                line = line.strip("\s\n")
                attrs = [float(attr) for attr in re.split("[,\s]+", line) if not attr == '']
                node_attrs.append(np.array(attrs))
    except IOError:
        print('No node attributes')
       
    label_has_zero = False
    filename_graphs=prefix + '_graph_labels.txt'
    graph_labels=[]
    with open(filename_graphs) as f:
        for line in f:
            line=line.strip("\n")
            val = int(line)
            if val == 0:
                label_has_zero = True
            graph_labels.append(val - 1)
    graph_labels = np.array(graph_labels)
    if label_has_zero:
        graph_labels += 1
    
    filename_adj=prefix + '_A.txt'
    adj_list={i:[] for i in range(1,len(graph_labels)+1)}    
    # index_graph={i:[] for i in range(1,len(graph_labels)+1)}
    num_edges = 0
    with open(filename_adj) as f:
        for line in f:
            line=line.strip("\n").split(",")
            e0,e1=(int(line[0].strip(" ")),int(line[1].strip(" ")))
            adj_list[graph_indic[e0]].append((e0,e1))
            # index_graph[graph_indic[e0]]+=[e0,e1]
            num_edges += 1
    # for k in index_graph.keys():
        # index_graph[k]=[u-1 for u in set(index_graph[k])]


    graphs=[]
    for i in range(1,1+len(adj_list)):
        # indexed from 1 here
        G=nx.from_edgelist(adj_list[i])
        graphs.append(G)
      
        # add features and labels
    for nodeid, nl in enumerate(node_labels):
        nodeid += 1
        graphs[graph_indic[nodeid]-1].add_node(nodeid)
        # graphs[graph_indic[nodeid]-1][nodeid]['label'] = nl

    for idx, G in enumerate(graphs):
        # no graph labels needed
        G.graph['label'] = graph_labels[idx]
        for u in G.nodes():
            if len(node_labels) > 0:
                G.node[u]['label'] = node_labels[u-1]
            if len(node_attrs) > 0:
                G.node[u]['feat'] = node_attrs[u-1]

        graphs[idx] = G

    # relabeling
    for idx, G in enumerate(graphs):
        mapping={}
        it=0
        if float(nx.__version__)<2.0:
            for n in G.nodes():
                mapping[n]=it
                it+=1
        else:
            for n in G.nodes:
                mapping[n]=it
                it+=1
            
        # indexed from 0
        G = nx.relabel_nodes(G, mapping)

        if max_nodes and G.number_of_nodes() > max_nodes:
            G = G.subgraph([i for i in range(0, max_nodes)])

        graphs[idx] = G

    # either max_nodes is 0 or None or 
    # return [g for g in graphs if g.number_of_nodes() <= 2000]
    return graphs

#def load_data(ds_name):
#    f = open(DATA_DIR + "/%s.graph"%ds_name, "r")
#    data = pickle.load(f)
#    graph_data = data["graph"]
#    labels = data["labels"]
#    labels  = np.array(labels, dtype = np.float)
#    print "Saving the label data into: "+ OUTPUT_DIR + "/%s_labels.mat"%ds_name
#    scipy.io.savemat(OUTPUT_DIR + "/%s_labels.mat"%ds_name, mdict={'label': labels})
#    return graph_data

def load_graph(graph):
    max_size = []
    for i in graph.nodes():
        max_size.append(i)
    for nidx in graph.nodes():
        for neighbor in graph.neighbors(nidx):
            max_size.append(neighbor)
    size = max(max_size)+1
    am = np.zeros((size, size))
    for nidx in graph.nodes():
        for neighbor in graph.neighbors(nidx):
            am[nidx][neighbor] = 1
    return am

def build_graphlet_corpus(ds_name, num_graphlets, samplesize):
    # if no graphlet is found in a graph, we will fall back to 0th graphlet of size k
    fallback_map = {1: 1, 2: 2, 3: 4, 4: 8, 5: 19, 6: 53, 7: 209, 8: 1253, 9: 13599}
    data = load_data(ds_name)
    len_data = len(data)
    canonical_map, weight_map = get_maps(num_graphlets)
    stat_arr = []
    vocabulary = set()
    corpus = []
    # randomly sample graphlets
    graph_map = {}
    graph_map2 = {}
    for gidx in range(len(data)):
        am = load_graph(data[gidx])
        m = len(am)
        #m = am.number_of_nodes()
        count_map = {}
        base_map = {}
        tmp_corpus = []
        if m >=num_graphlets:
            for j in range(samplesize):
                r =  np.random.permutation(range(m))
                for n in [num_graphlets]:
                    window = am[np.ix_(r[0:n],r[0:n])]
                    g_type = canonical_map[get_graphlet(window, n)]
                    graphlet_idx = g_type["idx"]
                    level = g_type["n"]
                    count_map[graphlet_idx] = count_map.get(graphlet_idx, 0) + 1
                    # for each node, pick a node it is connected to, and place a non-overlapping window
                    tmp_corpus.append(graphlet_idx)
                    vocabulary.add(graphlet_idx)
                    for node in r[0:n]:
                        # place a window to for each node in the original window
                        new_n_arr = r[n:][0:n-1] # select n-1 nodes
                        r2 = np.array(list(new_n_arr) + [node])
                        window2 = am[np.ix_(r2,r2)]
                        g_type2 = canonical_map[get_graphlet(window2, n)]
                        graphlet_idx2 = g_type2["idx"]
                        count_map[graphlet_idx2] = count_map.get(graphlet_idx2, 0) + 1
                        vocabulary.add(graphlet_idx2)
                        tmp_corpus.append(graphlet_idx2)
            corpus.append(tmp_corpus)
        else:
            count_map[fallback_map[num_graphlets]] = samplesize # fallback to 0th node at that level
        graph_map[gidx] = count_map
        #print "Graph: %s #nodes: %s  total samples: %s"%(gidx, len(data[gidx]), sum(graph_map[gidx].values()))

    #print "Total size of the corpus: ", len(corpus)
    prob_map = {gidx: {graphlet: count/float(sum(graphlets.values())) \
        for graphlet, count in graphlets.iteritems()} for gidx, graphlets in graph_map.iteritems()}
    num_graphs = len(prob_map)
    return corpus, vocabulary, prob_map, num_graphs

def get_graphlet(window, nsize):
    """
    This function takes the upper triangle of a nxn matrix and computes its canonical map
    """
    adj_mat = {idx: [i for i in list(np.where(edge)[0]) if i!=idx] for idx, edge in enumerate(window)}

    g = pynauty.Graph(number_of_vertices=nsize, directed=False, adjacency_dict = adj_mat)
    cert = pynauty.certificate(g)
    return cert

def get_maps(n):
    # canonical_map -> {canonical string id: {"graph", "idx", "n"}}
    file_counter = open("canonical_maps/canonical_map_n%s.p"%n, "rb")
    canonical_map = pickle.load(file_counter)
    file_counter.close()
    # weight map -> {parent id: {child1: weight1, ...}}
    file_counter = open("graphlet_counter_maps/graphlet_counter_nodebased_n%s.p"%n, "rb")
    weight_map = pickle.load(file_counter)
    file_counter.close()
    weight_map = {parent: {child: weight/float(sum(children.values())) for child, weight in children.iteritems()} for parent, children in weight_map.iteritems()}
    child_map = {}
    for parent, children in weight_map.iteritems():
        for k,v in children.iteritems():
            if k not in child_map:
                child_map[k] = {}
            child_map[k][parent] = v
    weight_map = child_map
    return canonical_map, weight_map

def adj_wrapper(g):
    am_ = g["al"]
    size =  max(np.shape(am_))
    am = np.zeros((size, size))
    for idx, i in enumerate(am_):
        for j in i:
            am[idx][j-1] = 1
    return am

def build_sp_corpus(ds_name, kernel_type=3):
    graph_data = load_data(ds_name)
    vocabulary = set()
    prob_map = {}
    graphs = {}
    corpus = []
    # compute MLE seperately for time-measuring purposes
    for gidx, graph in enumerate(graph_data):
        prob_map[gidx] = {}
        graphs[gidx] = []
        #label_map = [graph[nidx]["label"] for nidx in sorted(graph.keys())]
        label_map = [graph.node[nidx]["label"] for nidx in graph.nodes()]
        # construct networkX graph
        G = graph
        """
        G=nx.Graph()
        for nidx in graph:
            edges = graph[nidx]["neighbors"]
            if len(edges) == 0: # some nodes don't have any neighbors
                continue
            for edge in edges:
                G.add_edge(nidx, edge)
        """
        # get all pairs shortest paths
        all_shortest_paths = nx.all_pairs_shortest_path(G) # nx.floyd_warshall(G)
        # traverse all paths and subpaths
        tmp_corpus = []
        for source, sink_map in all_shortest_paths:#.iteritems():
            for sink, path in sink_map.iteritems():
                sp_length = len(path)-1
                #label = "_".join(map(str, sorted([label_map[source][0],label_map[sink][0]]))) + "_" + str(sp_length)
                label = "_".join(map(str, sorted([label_map[source],label_map[sink]]))) + "_" + str(sp_length)
                tmp_corpus.append(label)
                prob_map[gidx][label] = prob_map[gidx].get(label, 0) + 1
                graphs[gidx].append(label)
                vocabulary.add(label)
        corpus.append(tmp_corpus)
    if kernel_type != 3:
        # note that in MLE version, prob_map is not normalized.
        prob_map = {gidx: {path: count/float(sum(paths.values())) \
            for path, count in paths.iteritems()} for gidx, paths in prob_map.iteritems()}
    num_graphs = len(prob_map)
    return corpus, vocabulary, prob_map, num_graphs

def get_label(s):
    x = sorted([s[0], s[-1]])
    x = str(x[0][0]) + "_" + str(x[-1][0]) + "_" + str(len(s))
    return x

def build_wl_corpus(ds_name, max_h):
    graph_data = load_data(ds_name)
    labels = {}
    label_lookup = {}
    label_counter = 0
    vocabulary = set()
    num_graphs = len(graph_data)
    max_window = []
    wl_graph_map = {it: {gidx: defaultdict(lambda: 0) for gidx in range(num_graphs)} for it in range(-1, max_h)}
    sim_map = {}

    # initial labeling
    for gidx in range(num_graphs):
        #labels[gidx] = np.zeros(len(graph_data[gidx]), dtype = np.int32)
        labels[gidx] = np.zeros(graph_data[gidx].number_of_nodes(), dtype = np.int32)
        for node in range(len(graph_data[gidx])):
            label = graph_data[gidx].node[node]['label']
            if not label_lookup.has_key(label):
                label_lookup[label] = label_counter
                labels[gidx][node] = label_counter
                label_counter += 1
            else:
                labels[gidx][node] = label_lookup[label]
            wl_graph_map[-1][gidx][label_lookup[label]] = wl_graph_map[-1][gidx].get(label_lookup[label], 0) + 1
    compressed_labels = copy.deepcopy(labels)
    # WL iterations started
    for it in range(max_h):
        label_lookup = {}
        label_counter = 0
        for gidx in range(num_graphs):
            for node in graph_data[gidx].nodes():
                node_label = tuple([labels[gidx][node]])
                neighbors = list(graph_data[gidx].neighbors(node))
                if len(neighbors) > 0:
                    neighbors_label = tuple([labels[gidx][i] for i in neighbors])
                    node_label = tuple(tuple(node_label) + tuple(sorted(neighbors_label)))
                if not label_lookup.has_key(node_label):
                    label_lookup[node_label] = str(label_counter)
                    compressed_labels[gidx][node] = str(label_counter)
                    label_counter += 1
                else:
                    compressed_labels[gidx][node] = label_lookup[node_label]
                wl_graph_map[it][gidx][label_lookup[node_label]] = wl_graph_map[it][gidx].get(label_lookup[node_label], 0) + 1
        print "Number of compressed labels at iteration %s: %s"%(it, len(label_lookup))
        labels = copy.deepcopy(compressed_labels)

    # merge the following code into the loop above
    graphs = {}
    prob_map = {}
    corpus = []
    for it in range(-1, max_h):
        for gidx, label_map in wl_graph_map[it].iteritems():
            if gidx not in graphs:
                graphs[gidx] = []
                prob_map[gidx] = {}
            for label_, count in label_map.iteritems():
                label = str(it) + "+" + str(label_)
                for c in range(count):
                    graphs[gidx].append(label)
                vocabulary.add(label)
                prob_map[gidx][label] = count

    corpus = [graph for gidx, graph in graphs.iteritems()]
    vocabulary = sorted(vocabulary)
    return corpus, vocabulary, prob_map,  num_graphs

def l2_norm(vec):
    return  np.sqrt(np.dot(vec, vec))

if __name__ == "__main__":
    # location to save the results
    OUTPUT_DIR = "kernels/"
    # location of the datasets
    DATA_DIR = "datasets/"

    # hyperparameters
    num_dimensions = int(sys.argv[1]) # any integer > 0
    kernel_type = int(sys.argv[2]) # 1 (deep, l2) or 2 (deep, M), 3 (MLE)
    feature_type = int(sys.argv[3]) # 1 (graphlet), 2 (SP), 3 (WL)
    ds_name = sys.argv[4] # dataset name
    window_size = int(sys.argv[5]) # any integer > 0

    # graph kernel parameters
    max_h = 2 # for WL

    # word2vec parameters
    ngram_type = int(sys.argv[6]) # 1 (skip-gram), 0 (cbow)
    sampling_type = int(sys.argv[7]) # 1 (hierarchical sampling), 0 (negative sampling)
    graphlet_size = int(sys.argv[8]) # any integer > 0
    sample_size = int(sys.argv[9]) # any integer > 0
    print "Dataset: %s\n\nWord2vec Parameters:\nDimension: %s\nWindow size: %s\nNgram type: %s\nSampling type: %s\n\
            \n\nKernel-related Parameters:\nKernel type: %s\nFeature type: %s\nWL height: %s\nGraphlet size: %s\nSample size: %s\n"\
            %(ds_name, num_dimensions, window_size, ngram_type, sampling_type, kernel_type, feature_type, max_h, graphlet_size, sample_size)

    # STEP 1: Build corpus
    start = time.time()
    if feature_type == 1:
        # terms are graphlets
        corpus, vocabulary, prob_map, num_graphs = build_graphlet_corpus(ds_name, graphlet_size, sample_size)
    elif feature_type == 2:
        # terms are labeled shortest paths
        corpus, vocabulary, prob_map, num_graphs = build_sp_corpus(ds_name, kernel_type)
    elif feature_type == 3:
        # terms are labeled sub-trees
        corpus, vocabulary, prob_map, num_graphs = build_wl_corpus(ds_name, max_h)
    else:
        raise Exception("Unknown feature type!")
    end = time.time()
    vocabulary = list(sorted(vocabulary))
    print "Corpus construction total time: %g vocabulary size: %s"%(end-start, len(vocabulary))

    # STEP 2: learn hidden representations
    start = time.time()
    for i in range(len(corpus)):
        corpus[i] = list(map(str, corpus[i]))
    model = Word2Vec(corpus, size=num_dimensions, window=window_size, min_count=0, sg=ngram_type, hs=sampling_type)
    end = time.time()
    print "M matrix total time: %g"%(end-start)

    # STEP 3: compute the kernel
    K = np.zeros((num_graphs, num_graphs))
    if kernel_type == 1:
        # deep w/l2 norm
        P = np.zeros((num_graphs, len(vocabulary)))
        for i in range(num_graphs):
            for jdx, j in enumerate(vocabulary):
                P[i][jdx] = prob_map[i].get(j,0)
        M = np.zeros((len(vocabulary), len(vocabulary)))
        for idx,i in enumerate(vocabulary):
            M[idx][idx] = l2_norm(model[i])
        K = (P.dot(M)).dot(P.T)
    elif kernel_type == 2:
        P = np.zeros((num_graphs, len(vocabulary)))
        for i in range(num_graphs):
            for jdx, j in enumerate(vocabulary):
                P[i][jdx] = prob_map[i].get(j,0)
        M = np.zeros((len(vocabulary), len(vocabulary)))
        for i in range(len(vocabulary)):
            for j in range(len(vocabulary)):
                M[i][j] = np.dot(model[vocabulary[i]], model[vocabulary[j]])
        K = (P.dot(M)).dot(P.T)
    elif kernel_type == 3:
        #MLE
        P = np.zeros((num_graphs, len(vocabulary)))
        for i in range(num_graphs):
            for jdx, j in enumerate(vocabulary):
                P[i][jdx] = prob_map[i].get(j,0)
        K = P.dot(P.T)

    print(type(K))
    print(K.shape)
    # the following computed kernel can be directly fed to libsvm library
    # print "Saving the kernel to the following location: %s/deep_kernel_%s_k%s_d%s_f%s_w%s_ngram%s_sampling%s_gsize%s_samplesize%s.mat"%(OUTPUT_DIR, ds_name, kernel_type, num_dimensions, feature_type, window_size, ngram_type, sampling_type, graphlet_size, sample_size)
    # scipy.io.savemat("%s/deep_kernel_%s_k%s_d%s_f%s_w%s_ngram%s_sampling%s_gsize%s_samplesize%s.mat"%(OUTPUT_DIR, ds_name, kernel_type, num_dimensions, feature_type, window_size, ngram_type, sampling_type, graphlet_size, sample_size), mdict={'kernel': K})

    from sklearn.model_selection import GridSearchCV, StratifiedKFold
    from sklearn.svm import SVC, LinearSVC
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import accuracy_score


    graphs = load_data(ds_name)
    print('number of graphs', len(graphs))
    labels = [graph.graph['label'] for graph in graphs]
    
    x = np.array(graphs)
    y = np.array(labels)
    

    kf = StratifiedKFold(n_splits=10, random_state=None)
    kf.shuffle=True

    accuracies = []
    for train_index, test_index in kf.split(x, y):
        best_acc1 = 0
        x_train, x_test = x[train_index], x[test_index]
        K_train = K[train_index, :]
        K_train = K_train[:, train_index]
        #print(K_train.shape)

        y_train, y_test = y[train_index], y[test_index]
        K_test = K[test_index, :]
        K_test = K_test[:, train_index]
        #print(K_test.shape)

        # x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.1)
        search=True
        if search:
            params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
            classifier = GridSearchCV(SVC(kernel='precomputed'), params, cv=10, scoring='accuracy', verbose=0)
        else:
            classifier = SVC(C=10)
        classifier.fit(K_train, y_train)
        accuracies.append(accuracy_score(y_test, classifier.predict(K_test)))

    feature_type = int(sys.argv[3]) # 1 (graphlet), 2 (SP), 3 (WL)
    if feature_type == 1:
        tpe = 'graphlet'
    if feature_type == 2:
        tpe = 'sp'
    if feature_type == 3:
        tpe = 'wl'

    res =  np.mean(accuracies)
    with open('log', 'a+') as f:
        f.write('{},{}\n'.format(tpe, res))

