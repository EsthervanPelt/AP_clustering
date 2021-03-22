# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:06:52 2021

@author: 20183816
"""
# import modules
from data import *
from util import *
from kmeans import *

# assign datafiles
fileMetadata = "GDSC_metadata.csv"
fileRMAExpression = "GDSC_RNA_expression.csv"

# read in data
data, gene_names = read_data(fileMetadata, fileRMAExpression);

# choose value for k
k = 4;

# choose random centroids and assign datapoints
clusters = initial_clustering(k, data);
data, clusters = assign_datapoints(data, clusters, dist = 0)

S = [0]
# # iterate 
nr_iterations = 20

for i in range(nr_iterations):
    data, clusters = assign_datapoints(data, clusters, dist = 0)
    
    non_zero_clusters = {}
    for key in clusters:
        if clusters[key][1] != []: non_zero_clusters[key] = clusters[key] 
    
    clusters = compute_centroid(data, non_zero_clusters)
    
    S.append(silhouette_coeff(data, clusters, dist = 0))

    if S[-1] == S[-2]: break
 
print('Silhouette coefficient:', S[-1], '\n')

for key in clusters:
    
    NB = 0
    BRCA = 0
    KIRC = 0
    COAD_READ = 0
    
    for datapoint in clusters[key][1]:
        if data[datapoint][1] == 'NB': NB += 1
        if data[datapoint][1] == 'BRCA': BRCA += 1
        if data[datapoint][1] == 'KIRC': KIRC += 1
        if data[datapoint][1] == 'COAD/READ': COAD_READ += 1
        
    print('Cluster', key, '\n',
          'NB:       ', NB, '\n', 
          'BRCA:     ', BRCA, '\n',
          'KIRC:     ', KIRC, '\n',
          'COAD/READ:', COAD_READ, '\n')