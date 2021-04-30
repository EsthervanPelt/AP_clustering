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

def find_edge_threshold(data: dict, perc = 0.1, dist = 0):
    """
    Function that finds calculates the correlation matrix, and find a threshold such that a certain percentage of the values is above it.
    
    Parameters
    ----------
    data : dict
        key-value pairs with cell line indices as keys and instance values as second value.
    perc : float
        percentage of values in the correlation matrix that needs to be above the threshold
    dist : boolean value
        indicates whether the correlation, or the correlation distance should be used.
        0 = correlation, 1 = correlation distance. The default is 0.

    Returns
    -------
    threshold : float
        value above which the correlation must be 
    corr_matrix : list of lists
        filled with floats that contains the correlation between different instances in the data matrix
    datapoints : list 
        contains the gene indices corresponding to the rows/columns in the correlation matrix

    """
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
    
    # plt.figure()
    # plt.plot(c, fc)
    # plt.xlabel('c')
    # plt.ylabel('f(c)')
    # plt.title('fraction of nodes f(c) with correlation coefficient >= c')
    # plt.show()
    
    threshold = c[fc.index(list(filter(lambda i: i<=perc, fc))[0])]

    return threshold, corr_matrix, datapoints

def construct_graph(data: dict, threshold: float, perc: float, dist = 0):
    """
    Constructs a graph of the datapoints, in which all instances are nodes.
    Edges are added if the correlation between two nodes is above the threshold value.
    This threshold value, if not provided, is computed using the function find_edge_threshold.

    Parameters
    ----------
    data : dict
        key-value pairs with cell line indices as keys and instance values as second value.
    threshold : float
        value above which the correlation between data instances must be for nodes to be connected
    perc : float
        percentage of values in the correlation matrix that needs to be above the threshold
    dist : boolean value
        indicates whether the correlation, or the correlation distance should be used.
        0 = correlation, 1 = correlation distance. The default is 0.

    Returns
    -------
    G : multigraph
        graph of datapoints

    """
    G = nx.MultiGraph()
        
    if len(data)>1:
        new_threshold, corr_matrix, datapoints = find_edge_threshold(data, perc, dist)
        if threshold == None: threshold = new_threshold
        
        for key in data:
            G.add_node(key, label = data[key][1], gene_values = data[key][2], contains = []) #0 = cancer label, 1 = gene values, 2 = list of contained nodes
        if len(data) > 1:
            for row in range(len(corr_matrix)):
                for column in range(len(corr_matrix[row])):
                    if corr_matrix[row][column] >= threshold: 
                        G.add_edge(datapoints[row], datapoints[column])
                        
    elif len(data) == 1:
        for key in data:
            G.add_node(key, label = data[key][1], gene_values = data[key][2], contains = []) #0 = cancer label, 1 = gene values, 2 = list of contained nodes
        
    return G, threshold

def singleton(G):
    """
    Function that removes singeltons (i.e. nodes without edges) from a graph and bundles them into a new singleton graph.

    Parameters
    ----------
    G : multigraph
        Original graph with singletons

    Returns
    -------
    G : multigraph
        original graph, but now without singletons.
    S : graph
        set of singletons.

    """
    S = nx.Graph()
    
    to_remove = []
    
    for node in G.nodes():
        if G.degree[node] == 0:
            S.add_node(node)
            to_remove.append(node)
    
    for node in to_remove:
        G.remove_node(node)
        
    return G, S

def karger_min_cut(data: dict, G, threshold: float, dist = 0):
    """
    Function that performs a karger min cut of a graph. 
    This is done by contracting two random nodes, until the graph consists of two supernodes.
    Those supernodes then contain a number of contracted nodes. These two sets are used to constructing two new graphs.

    Parameters
    ----------
    data : dict
        key-value pairs with cell line indices as keys and instance values as second value.
    G : multigraph
        graph to perform the karger cut on, all nodes need to be connected to eachother
    threshold : float
        value above which the correlation between instances must be for nodes to be connected in the newly formed graphs
    dist : boolean value
        indicates whether the correlation, or the correlation distance should be used.
        0 = correlation, 1 = correlation distance. The default is 0.

    Returns
    -------
    mincut : integer
        number of edges that connect the two subsets
    listH : list
        list that contains the two newly constructed graphs, based on the two subsets from the karger cut
    list_data : list
        list that contains two dictionaries, containing the datasets for the two graphs 

    """

    while G.number_of_nodes() > 2:
        G = contract(G)
    
    for supernode in list(G.nodes): G.nodes[supernode]['contains'] += [supernode]
    mincut = len(G.edges(list(G.nodes())[0]))
    
    data1 = {key:value for key,value in data.items() if key in list(G.nodes().data('contains'))[0][1]}
    data2 = {key:value for key,value in data.items() if key in list(G.nodes().data('contains'))[1][1]}
    list_data = [data1, data2]
    
    listH = []
    if len(data1) > 0: H1, _ = construct_graph(data1, threshold, dist); listH.append(H1)
    if len(data2) > 0: H2, _ = construct_graph(data2, threshold, dist); listH.append(H2)
    if (len(data1) == 0) or len(data2) == 0: print("error: one node is empty")

    return mincut, listH, list_data
    
def contract(G):
    """
    Function that is used in the karger min cut algorithm, to contract two nodes.
    It adds the second node, neighbour one to the list of nodes that is contained in the first node, the main node.
    All edges that the neighbouring node was part of, are now connected to the main node.

    Parameters
    ----------
    G : multigraph
        graph containing nodes to be contracted

    Returns
    -------
    G : multigraph
        graph with one less node than the original graph, since two nodes are contracted

    """
    main_node = random.choice(list(G.nodes))
    neighbour = random.choice(list(G.adj[main_node]))
    
    for edge in G.edges(neighbour):
        if edge[0] == neighbour:
            if main_node != edge[1]:
                G.add_edge(main_node, edge[1])
        elif edge[1] == neighbour:
            if main_node != edge[0]:
                G.add_edge(main_node, edge[0])
    
    G.nodes[main_node]['contains'] += G.nodes[neighbour]['contains'] + [neighbour]
    
    G.remove_node(neighbour)
    return G

def DFSUtil(G, temp, v, visited, nodes):
    """
    Recursive function that creates a new graph from a set of connected nodes within the original graph.
    It marks the connected nodes from the old graph as visited.

    Parameters
    ----------
    G : multigraph
        graph from which the connected components should be retrieved.
    temp : multigraph
        new graph of the set of connected nodes
    v : integer
        index of node that is considered and added to the new graph
    visited : list of booleans
        has as length the number of nodes in the original nodes, contains true or false indicating whether each node is visited or not
    nodes : list
        list of nodes, the indices correspond to those of the 'visited' list

    Returns
    -------
    temp : multigraph
        new graph of the set of connected nodes

    """
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
    """
    Function that retrieves connected components in an undirected graph.

    Parameters
    ----------
    G : multigraph
        graph from which the connected components should be retrieved

    Returns
    -------
    cc : list
        list of multigraphs, which are the different connected components that were in G
        
    """
    visited = []
    cc = []
    for i in range(len(G.nodes)):
        visited.append(False)
      
    nodes = list(G.nodes)
    for v in range(len(nodes)):

        if visited[v] == False:
            temp = nx.MultiGraph()
            cc.append(DFSUtil(G, temp, v, visited, nodes))
    return cc

def adapted_hcs(data: dict, G, subgraphs: list, singles: list, threshold: float, mincut_trials: int, dist = 0):
    """
    Recursive function that transforms an initial (fully connected) (sub)graph to a highly connected subgraph.
    It first checks whether the input graph is already higher connected.
    If not, it performs several iterations of the karger cut algorithm until a minimal cut is obtained.
    This min cut set then forms two new subgraphs which are recursively put into this function to check whether they are highly connected too. 

    Parameters
    ----------
    data : dict
        key-value pairs with cell line indices as keys and instance values as second value.
    G : multigraph
        input graph, to be divided into highly connected subgraphs
    subgraphs : list
        will contains the highly connected subgraphs that have been formed
    singles : list
        if a newly formed subgraph consists of a single node, it is added to this list
    threshold : float
        value above which the correlation between instances must be for nodes to be connected in the newly formed graphs
    mincut_trials : int
        amount of times the karger cut algorithm should be run to find a minimal cut
    dist : dist : boolean value
        indicates whether the correlation, or the correlation distance should be used in the karger cut algorithm.
        0 = correlation, 1 = correlation distance. The default is 0.

    Returns
    -------
    subgraphs : TYPE
        contains the highly connected subgraphs that have been formed
    singles : TYPE
        if a newly formed subgraph consists of a single node, it is added to this list
    
    """
    highly_connected = True
    for node in G.nodes():
        if G.degree[node] < 0.5*G.number_of_nodes(): highly_connected = False;
    
    if (highly_connected==False) and (G.number_of_nodes() != 1):
        mincut, listH, list_data = karger_min_cut(data, G, threshold, dist)
            
        for i in range(mincut_trials-1):
            new_mincut, new_listH, new_list_data = karger_min_cut(data, G, threshold, dist)
            
            if new_mincut < mincut: 
                mincut = new_mincut; listH = new_listH; list_data = new_list_data
        
        subgraphs, singles = adapted_hcs(list_data[0], listH[0], subgraphs, singles, threshold, mincut_trials, dist)
        subgraphs, singles = adapted_hcs(list_data[1], listH[1], subgraphs, singles, threshold, mincut_trials, dist)
        
    elif (highly_connected == True) and (G.number_of_nodes() != 1): subgraphs.append(G)
    elif G.number_of_nodes == 1: singles.append(G)
    
    return subgraphs, singles