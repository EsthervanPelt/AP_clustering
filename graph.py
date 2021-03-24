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

def singleton(G):
    S = nx.Graph()
    
    to_remove = []
    
    for node in G.nodes():
        if G.degree[node] == 0:
            S.add_node(node)
            to_remove.append(node)
    
    for node in to_remove:
        G.remove_node(node)
        
    return G, S

def karger_min_cut(data: dict, G, cuts, dist = 0):
    i = 0
    while G.number_of_nodes() > 2:
        # print(i)
        # i += 1
        node = random.choice(list(G.nodes))
        while len(G.adj[node])<2:
            node = random.choice(list(G.nodes))
            # print(node)
        neighbour = random.choice(list(G.adj[node]))
        G = contract(G, node, neighbour)
        
    mincut = len(G.edges(list(G.nodes())[0]))
    
    cuts.append(mincut)
    
    data1 = {key:value for key,value in data.items() if key in list(G.nodes().data('contains'))[0][1]}
    data2 = {key:value for key,value in data.items() if key in list(G.nodes().data('contains'))[1][1]}
    
    H1 = construct_graph(data1, dist)
    H2 = construct_graph(data2, dist)
    
    return cuts, H1, H2, data1, data2
    
def contract(G, main_node, neighbour):
    for edge in G.edges(neighbour):
        if edge[0] == neighbour:
            if main_node != edge[1]:
                G.add_edge(main_node, edge[1])
        elif edge[1] == neighbour:
            if main_node != edge[0]:
                G.add_edge(main_node, edge[0])
    
    if G.nodes[neighbour]['contains'] != 0:
        G.nodes[main_node]['contains'] += G.nodes[neighbour]['contains']
    G.nodes[main_node]['contains'] += [neighbour]
    
    # print(main_node, neighbour)
    for i in list(G.nodes):
        if len(G.adj[i]) == 0: print(main_node, neighbour, i, 'warning', len(G.edges(neighbour)))
    # print(len(G.adj[main_node]), len(G.adj[neighbour]))
    
    G.remove_node(neighbour)
    return G

def adapted_hcs(data: dict, G, cuts: list, dist = 0):
    for node in G.nodes():
        if G.degree[node] < 0.5*G.number_of_nodes(): highly_connected = False; break
    
    if not(highly_connected) and G.number_of_nodes() != 1:
        C, H1, H2, data1, data2 = karger_min_cut(data, G, cuts, dist)
    
        p1 = adapted_hcs(data1, H1, C, dist)
        p2 = adapted_hcs(data2, H2, C, dist)
                
        G = nx.Graph()
        G.add_edge(p1, p2)
        print('hi')
    return G