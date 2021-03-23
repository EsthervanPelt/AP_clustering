# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:19:41 2021

@author: 20183816
"""
from util import half_correlation_matrix
import matplotlib.pyplot as plt
import networkx as nx
import random

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
    
    G = nx.MultiGraph()
    for key in data:
        G.add_node(key, label = data[key][1], gene_values = data[key][2], contains = []) #0 = cancer label, 1 = gene values, 2 = list of contained nodes
        
    for row in range(len(corr_matrix)):
        for column in range(len(corr_matrix[row])):
            if corr_matrix[row][column] >= cc: 
                G.add_edge(datapoints[row], datapoints[column])
                
    return G
 
def karger_min_cut(G, cuts):
    i = 0
    while G.number_of_nodes() > 2:
        print(i)
        i += 1
        node = random.choice(list(G.nodes))
        
        neighbour = random.choice(list(G.adj[node]))
        G = contract(G, node, neighbour)
        
    mincut = len(G.edges(list(G.nodes())[0]))
    
    cuts.append(mincut)
     
    H1
    
    return cuts#, H1, H2
    
def contract(G, main_node, neighbour):
    for edge in G.edges(neighbour):
        if main_node != edge[1]:
            G.add_edge(main_node, edge[1])
    
    
    if G.nodes[neighbour]['contains'] != 0:
        G.nodes[main_node]['contains'] += G.nodes[neighbour]['contains']
    G.nodes[main_node]['contains'] += [neighbour]
    
    G.remove_node(neighbour)
    
    return G

def adaptedHCS(G):
    if not(G.number_of_edges() > 0.5*G.number_of_nodes()) or (G.number_of_nodes() == 1):
        C, (H1, H2) = nx.minimum_cut(G, )
        
        p1 = adaptedHCS(H1)
        p2 = adaptedHCS(H2)
        
        G = nx.Graph()
        G.add_edge(p1, p2)
    
    return G