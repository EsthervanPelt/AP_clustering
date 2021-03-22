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
k = 5;

# perform initial clustering and compute centroids
data, clusters = initial_clustering(k, data);
clusters = compute_centroid(data, clusters);

# iterate 
nr_iterations = 10
for i in range(nr_iterations):
    data, clusters = assign_datapoints(data, clusters, dist = 0)
    # clusters = compute_centroid(data, clusters)
    # print(i)
    # S = silhouette_coeff(data, clusters, dist = 0)
    # print(S)