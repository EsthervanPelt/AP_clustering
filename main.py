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

# assign datafiles
fileMetadata = "GDSC_metadata.csv"
fileRMAExpression = "GDSC_RNA_expression.csv"

# read in data
data, gene_names = read_data(fileMetadata, fileRMAExpression);

# k-means algorithm
# k = 4
# max_iterations = 20
# dist = 0 # euclidean
# data, clusters, S = clustering(data, k, max_iterations, dist)

# construct graph with approximately 0.1*n(n-1)/2 edges
dist = 0
G = construct_graph(data, dist)

# make singleton set
G,S = singleton(G)
# retrieve sets of connected nodes
cc = connectedComponents(G)

# HCS algorithm
cuts = []
for subG in cc:
    subG = adapted_hcs(data, subG, cuts, dist)
