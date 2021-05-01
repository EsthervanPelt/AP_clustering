# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:06:52 2021

@author: 20183816
"""
# import modules
from data import read_data
from kmeans import clustering
from graph import construct_graph
from graph import adapted_hcs
from graph import singleton
from graph import connectedComponents
from util import euclidean

# assign datafiles
fileMetadata = "GDSC_metadata.csv"
fileRMAExpression = "GDSC_RNA_expression.csv"

# read in data
data, gene_names = read_data(fileMetadata, fileRMAExpression);

##############################################################################
#################### test k-means algorithm ##################################

# settings for k-means
max_iterations = 20; dist = 1 # squared euclidean distance
trials = 5

# open output file
outf_1 = open('results kmeans clustering.txt', 'w')
outf_1.write('Settings:'+' \n '+'max iterations: '+str(max_iterations)+' \n '+'number of trials: '+str(trials)+' \n '+'distance metric: '+str(dist)+'\n')

# repeat for different values of k
for k in range(2,10):
    S = 0
    # get different (initial) clusterings and store the best
    for i in range(trials):
        new_data, new_clusters, new_S = clustering(data, k, max_iterations, dist)
        if new_S > S:
            data = new_data; clusters = new_clusters; S = new_S
    
    # write results to the output file
    outf_1.write('Number of clusters: '+str(k)+' \t '+'Silhoette coefficient: '+str(S)+'\n')
    outf_1.write(' Division of data types per cluster:'+'\n')
    for key in clusters:
        outf_1.write('  Cluster '+str(key)+'\n'+ 
              '  NB:        '+str(clusters[key][2]['NB'])+'\n'+ 
              '  BRCA:      '+str(clusters[key][2]['BRCA'])+'\n'+
              '  KIRC:      '+str(clusters[key][2]['KIRC'])+'\n'+
              '  COAD/READ: '+str(clusters[key][2]['COAD/READ'])+'\n')

outf_1.close()

##############################################################################
################### test hcs algorithm #######################################

# construct graph with approximately 0.1*n(n-1)/2 edges
dist = 1; perc = 0.1
G, threshold = construct_graph(data, None, perc, dist)

# make singleton set
G,S = singleton(G)
# retrieve sets of connected nodes
cc = connectedComponents(G) 

# HCS algorithm
mincut_trials = 100
subgraphs = []; singles = []

# printing statements can be used to keep track of results
for i in range(len(cc)):
    # print("Connected component", i)
    G = cc[i]; subgraphs = []; singles = [] 
    # print(" # nodes before:", G.number_of_nodes())
    subgraphs, singles = adapted_hcs(data, G, subgraphs, singles, threshold, mincut_trials, dist)
    # print(" # nodes after:", G.number_of_nodes(), "\n # subgraphs:", len(subgraphs), "\n # singles:", len(singles), "\n")

##############################################################################
################### classification ###########################################

# leave-one-out algorithm for k-means clustering
max_iterations = 20; dist = 1; kclusters = 5
kmeans_trials = [1, 3, 5, 7, 10, 15, 20]

# open file to write results to
outf_2 = open('results classification.txt', 'w')
outf_2.write('Settings:'+'\n'+' max iterations: '+str(max_iterations)+'\n '+'number of clusters: '+str(kclusters)+'\n '+'distance metric: '+str(dist)+'\n')

# repeat for different amounts of k nearest neighbours
for kmeans in kmeans_trials:
    error_score = 0
    
    # pick one of the data points for the test set
    for test in data.keys():

        # assign the other data points to the training set
        training = {}
        for key in data.keys(): 
            if key != test: training[key] = data[key]
        
        # perform the k-means clustering
        training, clusters, S = clustering(training, kclusters, max_iterations, dist)
        
        # determine the k nearest neighbours
        distances = {}
        for key in training:
            test_dist = euclidean(data[test][2], training[key][2], dist)
            distances[key] = [test_dist, training[key][3]]
        nearest = sorted(distances.items(), key = lambda x: x[1][0])[0:kmeans]

        # determine the nearest cluster, based on the majority vote of the k nearest neighbours
        counts = {}
        for key in clusters:
            counts[key] = 0
            for i in nearest:
                if i[1][1] == key: counts[key] += 1
        
        test_cluster = sorted(counts.items(), key = lambda x: x[1], reverse = True)[0][0]

        # determine the nearest cluster label, depending on the majority vote of datapoints within the cluster
        cluster_type = sorted(clusters[test_cluster][2].items(), key = lambda x: x[1], reverse = True)[0][0]
        
        # if the cluster type doesn't correspond to the ground truth, add 1 to the error score
        if cluster_type != data[test][1]: error_score += 1; 

    # write results to ouput file
    outf_2.write('Number of k nearest neighbours: '+str(kmeans)+'\t'+'Error score: '+str(error_score)+'\n')
outf_2.close()
# leave-one-out algorithm for hcs
