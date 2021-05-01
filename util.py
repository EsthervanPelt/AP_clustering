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
    x : vector

    Returns
    -------
    mean of x

    """
    
    return sum(x)/len(x)

def euclidean(u: list, v: list, squared = 1):
    """
    Function that calculates and returns either the normal or the squared euclidean distance between two vectors.
    The two vectors must be of the same length.

    Parameters
    ----------
    u : list 
        vector 1
    v : list 
        vector 2
    squared : boolean value
        indicates whether the normal or squared euclidean distance should be returned.
        0 = normal, 1 = squared. The default is 1.

    Returns
    -------
    euc or squared_euc, depending on squared

    """
    # vectors must be of the same length
    if len(u) != len(v): print('Error: vectors not of the same length')
    
    elif len(u) == len(v):
        # calculate both distance metrics
        squared_euc = sum(map(lambda x, y: (x-y)**2, u, v))
        euc = squared_euc**(1/2)
        
        
        # return the right distanc
        if squared == 0: return euc
        elif squared == 1: return squared_euc
    
def manhattan(u: list, v: list):
    """
    Function that calculates and returns the manhattan distance between two vectors.
    The vectors must be of the same length.
    
    Parameters
    ----------
    u : list
        vector 1
    v : list
        vector 2

    Returns
    -------
    manh : calculated manhattan distance

    """
    # vectors must be of the same length
    if len(u) != len(v): print('Error: vectors not of the same length')
    
    elif len(u) == len(v):
        # calculate and return manhattan distance
        manh = sum(map(lambda x, y: abs(x-y), u, v))
        return manh

def correlation(u: list, v: list, dist = 1):
    """
    Function that calculates the correlation between two vectors. 
    It returns either the correlation, or the correlation distance (which is equal to 1 - |correlation|).

    Parameters
    ----------
    u : list 
        vector 1
    v : list 
        vector 2
    dist : boolean value
        indicates whether the correlation, or the correlation distance should be returned.
        0 = correlation, 1 = correlation distance. The default is 1.

    Returns
    -------
    corr or d_corr, depending on dist

    """
    # vectors must be of the same length
    if len(u) != len(v): print('Error: vectors not of the same length')
    
    elif len(u) == len(v):
        # calculate mean of both vectors
        mean_u = mean(u)
        mean_v = mean(v)
        
        # calculate nominator and denominator of correlation expression
        nom = sum(map(lambda x,y: (x-mean_u)*(y-mean_v), u, v))
        den = (sum([(x-mean_u)**2 for x in u])*sum([(y-mean_v)**2 for y in v]))**(1/2)
        
        # calculate both correlation and the correlation distance
        corr = nom/den 
        d_corr = 1 - abs(corr)
        
        # return either correlation or correlation distance
        if dist == 1: return d_corr
        elif dist == 0: return corr

def half_correlation_matrix(data: dict, dist = 0): 
    """
    Function that fills a correlation matrix in an efficient way such that only half of the values need to be computed

    Parameters
    ----------
    data : dict
        key-value pairs with data points as keys and values at the second place in the values tuple.
    dist : boolean value
        indicates whether the correlation, or the correlation distance should be used.
        0 = correlation, 1 = correlation distance. The default is 0.

    Returns
    -------
    corr_matrix : list of lists
        filled with floats that contains the correlation between different instances in the data matrix
    datapoints : list 
        contains the cell line indices corresponding to the rows/columns in the correlation matrix

    """
    # get a list of data points
    datapoints = list(data.keys())
    
    # create an empty list for the correlation matrix
    corr_matrix = []
    
    # loop over all instances
    for i in range(len(datapoints)):
        x = datapoints[i]
        
        row = []
        
        # loop over all datapoints up until the datapoint you're at
        for j in range (i):
            y = datapoints[j]
            # calculate correlation between the two datapoints
            corr = correlation(data[x][2], data[y][2], dist)
            # fill the current row
            row.append(corr)
        
        # add the row to the matrix
        corr_matrix.append(row)
    
    # return the correlation matrices and the list of data points
    return corr_matrix, datapoints

def mean_cluster_weight(x: int, data: dict, cluster: tuple, dist: 0): #0 = euclidean, 1 = squared euclidean
    """
    Function that calculates either the normal or the squared euclidean distance from one datapoint to all datapoints within a certain cluster.
    These distances are summed and then divided by the length of the cluster, to calculate the mean cluster weight.

    Parameters
    ----------
    x : int
        selected data point
    data : dict
        key-value pairs with data points as keys and values at the second place in the values tuple.
    cluster : tuple
        tuple that contains at the second place a list of data points that make up the cluster
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    MC : the mean cluster weight
        
    """
    # get list of data points within the cluster
    cluster_points = cluster[1]
    # get  values belonging to the selected data point
    x_values = data[x][2]
    
    sum_dist = 0
    # loop over all data points within the cluster
    for y in cluster_points:
        # get values of data point within the cluster
        y_values = data[y][2]
        
        # calculate distance, based on chosen distance metric
        if dist == 0: 
            sum_dist += euclidean(x_values, y_values, squared = 0)
        elif dist == 1:
            sum_dist += euclidean(x_values, y_values, squared = 1)
    
    # calculate the mean cluster weight
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
        selected data point
    data : dict
        key-value pairs with data points as keys and values at the second place in the values tuple.
    clusters : dict
        dictionary with as keys the cluster IDs and as values the lists of data points within the cluster
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    NC : the MC of the nearest cluster to the datapoint.
    
    """
    NC = None
    
    # loop over all clusters
    for cluster in clusters.keys():
        
        # check whether the selected data point is not within this cluster
        if cluster != data[x][3]:
            # calculate the mean cluster weight
            MC = mean_cluster_weight(x, data, clusters[cluster], dist)
            
            # if this is the first cluster weight calculated, or this is the smallest so far, assign this as the neighbouring cluster weight
            if type(MC) != None:
                if (NC == None) or (NC > MC): NC = MC
    
    return NC

def silhouette_score(x: int, data: dict, clusters: dict, dist: 0): #0 = euclidean, 1 = squared euclidean
    """
    Function that calculates the Silhouette score of a certain datapoint.
    
    Parameters
    ----------
    x : int
        selected data point 
    data : dict
        key-value pairs with data points as keys and values at the second place in the values tuple.
    clusters : dict
        dictionary with as keys the cluster IDs and as values the lists of data points within the cluster
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    s : the silhouette score
    
    """
    # get the cluster to which this data point belongs
    cluster_x = data[x][3]
            
    # get the mean cluster weight
    W = mean_cluster_weight(x, data, clusters[cluster_x], dist)
    # get the neighbouring cluster weight
    NC = neighbouring_cluster(x, data, clusters, dist) 
    
    # calculate the silhouette score
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
        key-value pairs with data points as keys and values at the second place in the values tuple.
    clusters : dict
        dictionary with as keys the cluster IDs and as values the lists of data points within the cluster
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    S : Silhouette coefficient
    
    """    
    # calculate the silhouette scores of all datapoints within the set
    sum_s = sum(map(lambda x: silhouette_score(x, data, clusters, dist), data.keys()))
    # calculate the silhouette coefficient
    S = sum_s / len(data.keys())
    
    return S

    
    