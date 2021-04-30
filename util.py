# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 09:51:15 2021

@author: 20183816
"""
def mean(x: list):
    """
    Function that calculates the mean of a list of values.
    
    Parameters
    ----------
    x : list of values

    Returns
    -------
    mean of x

    """
    
    return sum(x)/len(x)

def euclidean(u: list, v: list, squared = 1):
    """
    Function that calculates and returns either the normal or the squared euclidean distance between two lists of values.

    Parameters
    ----------
    u : list of values
    v : list of values
        u and v must be of the same length
    squared : boolean value
        indicates whether the normal or squared euclidean distance should be returned.
        0 = normal, 1 = squared. The default is 1.

    Returns
    -------
    euc or squared_euc, depending on squared

    """
    if len(u) != len(v): print('Error: vectors not of the same length')
    
    elif len(u) == len(v):
        squared_euc = sum(map(lambda x, y: (x-y)**2, u, v))
        euc = squared_euc**(1/2)
        
        if squared == 0: return euc
        elif squared == 1: return squared_euc
    
def manhattan(u: list, v: list):
    """
    Function that calculates and returns the manhattan distance between two lists of values.

    Parameters
    ----------
    u : list of values
    v : list of values

    Returns
    -------
    manh : calculated manhattan distance

    """
    if len(u) != len(v): print('Error: vectors not of the same length')
    
    elif len(u) == len(v):
        manh = sum(map(lambda x, y: abs(x-y), u, v))
        return manh

def correlation(u: list, v: list, dist = 1):
    """
    Function that calculates the correlation between two lists of values. 
    It returns either the correlation, or the correlation distance (which is equal to 1 - |correlation|).

    Parameters
    ----------
    u : list of values
    v : list of values
    dist : boolean value
        indicates whether the correlation, or the correlation distance should be returned.
        0 = correlation, 1 = correlation distance. The default is 1.

    Returns
    -------
    corr or d_corr, depending on dist

    """
    if len(u) != len(v): print('Error: vectors not of the same length')
    
    elif len(u) == len(v):
        mean_u = mean(u)
        mean_v = mean(v)
        
        nom = sum(map(lambda x,y: (x-mean_u)*(y-mean_v), u, v))
        den = (sum([(x-mean_u)**2 for x in u])*sum([(y-mean_v)**2 for y in v]))**(1/2)
        
        corr = nom/den 
        d_corr = 1 - abs(corr)
        
        if dist == 1: return d_corr
        elif dist == 0: return corr

def half_correlation_matrix(data: dict, dist = 0): 
    """
    Function that fills a correlation matrix in an efficient way such that only half of the values need to be computed

    Parameters
    ----------
    data : dict
        key-value pairs with cell line indices as keys and instance values as second value.
    dist : boolean value
        indicates whether the correlation, or the correlation distance should be used.
        0 = correlation, 1 = correlation distance. The default is 0.

    Returns
    -------
    corr_matrix : list of lists
        filled with floats that contains the correlation between different instances in the data matrix
    datapoints : list 
        contains the gene indices corresponding to the rows/columns in the correlation matrix

    """
    datapoints = list(data.keys())
    corr_matrix = []
    for i in range(len(datapoints)):
        x = datapoints[i]
        
        row = []
        
        for j in range (i):
            y = datapoints[j]
            corr = correlation(data[x][2], data[y][2], dist)
            row.append(corr)
        
        corr_matrix.append(row)
    
    return corr_matrix, datapoints

def mean_cluster_weight(x: int, data: dict, cluster: tuple, dist: 0): #0 = euclidean, 1 = squared euclidean
    """
    Function that calculates either the normal or the squared euclidean distance from one datapoint to all datapoints within a certain cluster.
    These distances are summed and then divided by the length of the cluster, to calculate the mean cluster weight.

    Parameters
    ----------
    x : int
        gene index of one datapoint 
    data : dict
        key-value pairs with cell line indices as keys and instance values as second value.
    cluster : tuple
        set of gene indices that together form one tuple
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    MC : the mean cluster weight
        
    """
    cluster_points = cluster[1]
    x_values = data[x][2]
    
    sum_dist = 0
    for y in cluster_points:
        y_values = data[y][2]
        
        if dist == 0: 
            sum_dist += euclidean(x_values, y_values, squared = 0)
        elif dist == 1:
            sum_dist += euclidean(x_values, y_values, squared = 1)
    
    MC = sum_dist / len(cluster_points)
    
    return MC

def neighbouring_cluster(x: int, data: dict, clusters: dict, dist = 0): #0 = euclidean, 1 = squared euclidean
    """
    Function that calculates the mean cluster weight of one datapoint with respect to all clusters.
    Then the cluster with the smallest MC is the neighbouring cluster, and this value is returned. 
    Either the normal, or the squared euclidean distance can be used.

    Parameters
    ----------
    x : int
        gene index of one datapoint 
    data : dict
        key-value pairs with cell line indices as keys and instance values as second value.
    clusters : dict
        dictionaries that contains the different sets of cluster of gene indices
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    NC : the MC of the nearest cluster to the datapoint.
    
    """
    NC = None
    for cluster in clusters.keys():
        if cluster != data[x][3]:
            MC = mean_cluster_weight(x, data, clusters[cluster], dist)
            
            if type(MC) != None:
                if (NC == None) or (NC > MC): NC = MC
    
    return NC

def silhouette_score(x: int, data: dict, clusters: dict, dist: 0): #0 = euclidean, 1 = squared euclidean
    """
    Function that calculates the Silhouette score of a certain datapoint.
    
    Parameters
    ----------
    x : int
        gene index of one datapoint 
    data : dict
        key-value pairs with cell line indices as keys and instance values as second value.
    clusters : dict
        dictionaries that contains the different sets of cluster of gene indices
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    s : the silhouette score
    
    """
    cluster_x = data[x][3]
            
    W = mean_cluster_weight(x, data, clusters[cluster_x], dist)
    NC = neighbouring_cluster(x, data, clusters, dist) 

    if W < NC: s = 1-W/NC
    elif W == NC: s = 0
    elif W > NC: s = NC/W-1
    
    return s

def silhouette_coeff(data: dict, clusters: dict, dist = 0): #0 = euclidean, 1 = squared euclidean
    """
    Function that calculates the Silhouette coefficient of a dataset, by calculating the mean silhouette score of all datapoints within the set. 

    Parameters
    ----------
    data : dict
        key-value pairs with cell line indices as keys and instance values as second value.
    clusters : dict
        dictionaries that contains the different sets of cluster of gene indices
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    S : Silhouette coefficient
    
    """    
    sum_s = sum(map(lambda x: silhouette_score(x, data, clusters, dist), data.keys()))
    S = sum_s / len(data.keys())
    
    return S

    
    