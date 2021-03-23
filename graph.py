# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:19:41 2021

@author: 20183816
"""
from util import half_correlation_matrix
import matplotlib.pyplot as plt
import networkx as nx
import random

####### Read in data ###########
from data import read_data
fileMetadata = "GDSC_metadata.csv"
fileRMAExpression = "GDSC_RNA_expression.csv"
data, gene_names = read_data(fileMetadata, fileRMAExpression);
################################

def find_edge_threshold(data: dict, dist = 0):
    corr_matrix, datapoints = half_correlation_matrix(data,dist=0)
    
    fc = []
    y = 200
    c = [x/(y-1) for x in range(y)]
    for cc in c:
        f = 0
        n = len(corr_matrix)
        for row in corr_matrix:
            for i in row: 
                if i >= cc: f +=1
                
        f = f/(n*(n-1)/2)
        
        fc.append(f)
    
    plt.figure()
    plt.plot(c, fc)
    plt.xlabel('c')
    plt.ylabel('f(c)')
    plt.title('fraction of nodes f(c) with correlation coefficient >= c')
    plt.show()
    
    cc = c[fc.index(list(filter(lambda i: i<=0.1, fc))[0])]

    return cc, corr_matrix, datapoints

def construct_graph(data: dict, dist = 0):
    cc, corr_matrix, datapoints = find_edge_threshold(data, dist)
    
    G = nx.Graph()
    for key in data:
        G.add_node(key, label = data[key][1], gene_values = data[key][2]) #0 = cancer label, 1 = gene values
    
    for row in range(len(corr_matrix)):
        for column in range(len(corr_matrix[row])):
            if corr_matrix[row][column] >= cc: 
                G.add_edge(datapoints[row], datapoints[column])
            
    return G
 
cuts = []
def kargerMinCut(G):
    while len(G) > 2:
        v = random.choice(list(G.keys()))
        w = random.choice(G[v])
        contract(G, v, w)
    mincut = len(G[list(G.keys())[0]])
    cuts.append(mincut)
    
def contract(graph, v, w):
    for node in graph[w]:  # merge the nodes from w to v
         if node != v:  # we dont want to add self-loops
             graph[v].append(node)
         graph[node].remove(w)  # delete the edges to the absorbed 
         if node != v:
              graph[node].append(v)
    del graph[w]  # delete the absorbed vertex 'w'

# def adaptedHCS(G = nx.Graph()):
#     if highly_connected(G):
#         return G
#     else:
#         H1, H2, C = KargerCut(G)
#         p1 = adaptedHCS(H1)
#         p2 = adaptedHCS(H2)