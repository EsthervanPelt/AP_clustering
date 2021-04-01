# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:19:41 2021

@author: 20183816
"""
from util import half_correlation_matrix
import matplotlib.pyplot as plt
import networkx as nx
import random
import matplotlib.pyplot

def find_edge_threshold(data: dict, dist = 0):
    corr_matrix, datapoints = half_correlation_matrix(data, dist=0)
    
    fc = []
    y = 200
    c = [x/(y-1) for x in range(y)]
    n = len(corr_matrix)
    for cc in c:
        f = 0
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
    
    threshold = c[fc.index(list(filter(lambda i: i<=0.1, fc))[0])]

    return threshold, corr_matrix, datapoints

def construct_graph(data: dict, threshold: float, dist = 0):
    new_threshold, corr_matrix, datapoints = find_edge_threshold(data, dist)
    if threshold == None: threshold = new_threshold
    
    G = nx.MultiGraph()
    for key in data:
        G.add_node(key, label = data[key][1], gene_values = data[key][2], contains = []) #0 = cancer label, 1 = gene values, 2 = list of contained nodes
        
    for row in range(len(corr_matrix)):
        for column in range(len(corr_matrix[row])):
            if corr_matrix[row][column] >= threshold: 
                G.add_edge(datapoints[row], datapoints[column])
    
    return G, threshold

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

def karger_min_cut(data: dict, G, cuts: list, threshold: float, dist = 0):

    while G.number_of_nodes() > 2:
        node = random.choice(list(G.nodes))
        
        check = True
        for i in G.nodes:
            if len(G.adj[i]) > 1: check = False
        if check==True: print('error'); break    
        
        while len(G.adj[node])<1:
            node = random.choice(list(G.nodes))
        neighbour = random.choice(list(G.adj[node]))
        G = contract(G, node, neighbour)
    print(G.nodes)
    print(G.edges)
    mincut = len(G.edges(list(G.nodes())[0]))
    
    cuts.append(mincut)
    
    data1 = {key:value for key,value in data.items() if key in list(G.nodes().data('contains'))[0][1]}
    data2 = {key:value for key,value in data.items() if key in list(G.nodes().data('contains'))[1][1]}
    
    listH = []
    if len(data1) > 1: listH.append(construct_graph(data1, threshold, dist))
    if len(data2) > 1: listH.append(construct_graph(data2, threshold, dist))
    
    return cuts, listH, data1, data2
    
def contract(G, main_node, neighbour):
    for edge in G.edges(neighbour):
        if edge[0] == neighbour:
            if main_node != edge[1]:
                G.add_edge(main_node, edge[1])
        elif edge[1] == neighbour:
            if main_node != edge[0]:
                G.add_edge(main_node, edge[0])
    
    G.nodes[main_node]['contains'] += G.nodes[neighbour]['contains'] + [neighbour]
    
    # print(main_node, neighbour)
    for i in list(G.nodes):
        if len(G.adj[i]) == 0: 
            print(main_node, neighbour, i, 'warning', len(G.edges(neighbour)))
    
    G.remove_node(neighbour)
    return G

def DFSUtil(G, temp, v, visited, nodes):
    node = nodes[v]
    
    # Mark the current vertex as visited
    visited[v] = True
 
    # Store the vertex to list

    attribute = G.nodes[node]
    temp.add_node(node, label = attribute['label'], gene_values = attribute['gene_values'], contains = attribute['contains'])
 
    # Repeat for all vertices adjacent to this vertex v
    
    for i in G.adj[node]:
        temp.add_edge(node, i)
        if visited[nodes.index(i)] == False:

            # Update the list
            temp = DFSUtil(G, temp, nodes.index(i), visited, nodes)
    return temp
 
# Method to retrieve connected components in an undirected graph
def connectedComponents(G):
    visited = []
    cc = []
    for i in range(len(G.nodes)):
        visited.append(False)
        
    nodes = list(G.nodes)
    for v in range(len(nodes)):

        if visited[v] == False:
            temp =nx.MultiGraph()
            cc.append(DFSUtil(G, temp, v, visited, nodes))
    return cc

def adapted_hcs(data: dict, G, cuts: list, threshold: float, dist = 0):
    highly_connected = True
    for node in G.nodes():
        if G.degree[node] < 0.5*G.number_of_nodes(): highly_connected = False;

    if (highly_connected==False) and G.number_of_nodes() != 1:
        cuts, listH, data1, data2 = karger_min_cut(data, G, cuts, threshold, dist)
        
        listp = []
        for H in listH:
            listp.append(adapted_hcs(data1, H, cuts, threshold, dist))
                
        G = nx.MultiGraph()
        if len(listp) > 1: G.add_edge(listp[0], listp[1])
        elif len(listp) == 1: G.add_node(listp[0])
        print('hi')
        
    return G, cuts