# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:19:41 2021

@author: 20183816
"""
from util import half_correlation_matrix
import matplotlib.pyplot as plt
import networkx as nx
import random

def find_edge_threshold(data: dict, perc = 0.1, dist = 0, plot = 0):
    """
    Function that finds calculates the correlation matrix, and find a threshold such that a certain percentage of the values is above it.
    
    Parameters
    ----------
    data : dict
        key-value pairs with data points as keys and instance values as second value.
    perc : float
        percentage of values in the correlation matrix that needs to be above the threshold
    dist : boolean value
        indicates whether the correlation, or the correlation distance should be used.
        0 = correlation, 1 = correlation distance. The default is 0.
    plot : boolean value
        determines whether a plot of f(c) should be made, with f the fraction of nodes above c
        0 = no plot, 1 = plot. The default is 0.
        
    Returns
    -------
    threshold : float
        value above which the correlation must be 
    corr_matrix : list of lists
        filled with floats that contains the correlation between different instances in the data matrix
    datapoints : list 
        contains the gene indices corresponding to the rows/columns in the correlation matrix

    """
    # calculate the correlation matrix
    corr_matrix, datapoints = half_correlation_matrix(data, dist=0)
    
    # create empty list to store fractions in
    fc = []
    # variant of linspace function, to create array of correlation coefficients
    y = 200; cc = [x/(y-1) for x in range(y)]
    
    # determine number of data points
    n = len(corr_matrix)
    
    # iterate over array of correlation coefficients
    for c in cc:
        f = 0
        # iterate over all points in the correlation matrix
        for row_index in range(n):
            for column_index in range(row_index): 
                i = corr_matrix[row_index][column_index]
                # count number of correlation coefficients above c
                if i >= c: f +=1
        
        # calculate fraction, account for the fact that the matrix is symmetric
        f = f/(n*(n-1)/2)
        # add fraction to list
        fc.append(f)
    
    # make a plot of the fraction above a certain correlation coefficient, against the coefficients
    if plot == 1:
        plt.figure()
        plt.plot(cc, fc)
        plt.xlabel('c')
        plt.ylabel('f(c)')
        plt.title('fraction of nodes f(c) with correlation coefficient >= c')
        plt.show()
    
    # determine the correlation threshold to obtain a given fraction
    threshold = cc[fc.index(list(filter(lambda i: i<=perc, fc))[0])]

    return threshold, corr_matrix, datapoints

def construct_graph(data: dict, threshold: float, perc: float, dist = 0):
    """
    Constructs a graph of the datapoints, in which all instances are nodes.
    Edges are added if the correlation between two nodes is above the threshold value.
    This threshold value, if not provided, is computed using the function find_edge_threshold.

    Parameters
    ----------
    data : dict
        key-value pairs with data points as keys and instance values as second value.
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
    # create empty graph, use multigraph to allow for multiple edges between nodes
    G = nx.MultiGraph()
    
    # if there are multiple data points
    if len(data)>1:
        # calculate threshold, correlation matrix and get list of datapoints
        new_threshold, corr_matrix, datapoints = find_edge_threshold(data, perc, dist, plot = 0)
        # if a threshold was not provided, use the newly calculated threshold
        if threshold == None: threshold = new_threshold
        
        # loop over all nodes and add to graph
        for key in data:
            G.add_node(key, label = data[key][1], gene_values = data[key][2], contains = []) #0 = cancer label, 1 = gene values, 2 = list of contained nodes
        
        # loop over half the correlation matrix (all combinations of data points)
        for row in range(len(corr_matrix)):
            for column in range(len(corr_matrix[row])):
                # if the correlation coefficient is above the threshold, add an edge between the two datapoins
                if corr_matrix[row][column] >= threshold: 
                    G.add_edge(datapoints[row], datapoints[column])
                        
    # if there is only one data point, create a graph with a single node
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
    # create empty graph
    S = nx.Graph()
    
    to_remove = []
    
    # check wheteher a node has neighbours
    for node in G.nodes():
        if G.degree[node] == 0:
            # if not, add the node to the singleton graph, and to the list of nodes to remove
            S.add_node(node)
            to_remove.append(node)
            
    # remove singletons from G
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
    # randomly contract nodes, until there are only two nodes left in G
    while G.number_of_nodes() > 2:
        G = contract(G)
    
    # add remaining nodes to their own lists of nodes contained
    for supernode in list(G.nodes): G.nodes[supernode]['contains'] += [supernode]
    # determine the quality of the cut
    mincut = len(G.edges(list(G.nodes())[0]))
    
    # split G into two new sets, containing all contracted nodes + the supernode of that set, make new data dictionaries
    data1 = {key:value for key,value in data.items() if key in list(G.nodes().data('contains'))[0][1]}
    data2 = {key:value for key,value in data.items() if key in list(G.nodes().data('contains'))[1][1]}
    # store new data dicts in one list
    list_data = [data1, data2]
    
    # create lists for new subgraphs to be stored in
    listH = []
    # construct new graphs from the new data dicts, using the threshold that has been calculated in constructing the original graph
    if len(data1) > 0: H1, _ = construct_graph(data1, threshold, dist); listH.append(H1)
    if len(data2) > 0: H2, _ = construct_graph(data2, threshold, dist); listH.append(H2)
    # print an error message if one set is empty
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
    # randomly choose two nodes from G
    main_node = random.choice(list(G.nodes))
    neighbour = random.choice(list(G.adj[main_node]))
    
    # loop over all edges of the neighbouring node
    for edge in G.edges(neighbour):
        # check whether it is not an edge between the two chosen nodes to avoid self loops
        if edge[0] == neighbour:
            if main_node != edge[1]:
                # if not, add the edge to the main node
                G.add_edge(main_node, edge[1])
        # check whether it is not an edge between the two chosen nodes to avoid self loops
        elif edge[1] == neighbour:
            if main_node != edge[0]:
                # if not, add the edge to the main node
                G.add_edge(main_node, edge[0])
    
    # store contracted node in list of contracted nodes of the main node
    G.nodes[main_node]['contains'] += G.nodes[neighbour]['contains'] + [neighbour]
    # remove contracted node
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
    # create empty list for visited nodes
    visited = []
    # create empty list for set of connected nodes
    cc = []
    
    # loop over nodes and mark them as unvisited
    for i in range(len(G.nodes)):
        visited.append(False)
      
    # loop over list of nodes
    nodes = list(G.nodes)
    for v in range(len(nodes)):

        # if not yet visited, check whether it is connected and make a subgraph of the connected set
        if visited[v] == False:
            temp = nx.MultiGraph() # create empty subgraph
            cc.append(DFSUtil(G, temp, v, visited, nodes)) # fill subgraph
            
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
    # initially assume the graph is highly connected
    highly_connected = True
    # loop over all nodes
    for node in G.nodes():
        # change the highly connected status, if the criterion is met
        if G.degree[node] < 0.5*G.number_of_nodes(): highly_connected = False;
    
    # if the graph consists of multiple nodes and is not highly connected
    if (highly_connected==False) and (G.number_of_nodes() != 1):
        # perform an initial karger cut
        mincut, listH, list_data = karger_min_cut(data, G, threshold, dist)
            
        # perform the rest of the karger cut
        for i in range(mincut_trials-1):
            new_mincut, new_listH, new_list_data = karger_min_cut(data, G, threshold, dist)
            
            # after each cut, check whether the cut is better then before and if so store the new cut 
            if new_mincut < mincut: 
                mincut = new_mincut; listH = new_listH; list_data = new_list_data
        
        # recursively perform this algorithm on the resulting cut set, to ensure that all subgraphs are highly connected
        subgraphs, singles = adapted_hcs(list_data[0], listH[0], subgraphs, singles, threshold, mincut_trials, dist)
        subgraphs, singles = adapted_hcs(list_data[1], listH[1], subgraphs, singles, threshold, mincut_trials, dist)
    
    # if the graph is highly connected, store it in the list of subgraphs
    elif (highly_connected == True) and (G.number_of_nodes() != 1): subgraphs.append(G)
    # if the graph consists of a single node, store it in the list of singles
    elif G.number_of_nodes == 1: singles.append(G)
    
    return subgraphs, singles