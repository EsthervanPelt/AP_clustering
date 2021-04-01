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
import networkx as nx
import matplotlib.pyplot as plt
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
G, threshold = construct_graph(data, None, dist)

# make singleton set
G,S = singleton(G)
# retrieve sets of connected nodes
cc = connectedComponents(G) #you can also use nx.connected_components()

# HCS algorithm
cuts = []
subG0 = cc[0]
plt.figure()
nx.draw(subG0)
subG1 = cc[1]
plt.figure()
nx.draw(subG1)
subG2 = cc[2]
plt.figure()
nx.draw(subG2)
subG, cuts = adapted_hcs(data, subG0, cuts, threshold, dist)
i = 0
# for subG in cc:
#     i+=1
#     print(i)
#     subG, cuts = adapted_hcs(data, subG, cuts, dist)
