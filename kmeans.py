# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 12:15:42 2021

@author: 20183816
"""
import random as rd
from util import euclidean
from util import silhouette_coeff

def initial_clustering(k: int, data: dict):
    """
    Function that performs a random initial clustering of the data points

    Parameters
    ----------
    k : int
        number of clusters to be generated
    data : dict
        key-value pairs with data points as keys and cluster ID at the third place in the values tuple (this now contains Nones as placeholders)

    Returns
    -------
    clusters : dict
        contains the cluster IDs as keys, and a list of cell line indices as value

    """
    # create an empty dictionary for the clusters to be stored in
    clusters = {}
    
    # create k clusters
    for cluster_ID in range(k):
        # choose a random data point 
        random_point = rd.choice(list(data.keys()))
        # use this chosen data point as initial centroid
        clusters[cluster_ID] = [data[random_point][2], [], {}] #0 = centroid location, 1 = list of data points, 2 = stats (data types)
        
    return clusters

def compute_centroid(data: dict, clusters: dict):
    """
    Function that computes the cluster centroid for each of the clusters. It uses the list of data points contained in the clusters dictionairy,
    combined with the values of those data points that are contained within the data dictionairy. 

    Parameters
    ----------
    data : dict
        key-value pairs with data points as keys and values at the second place in the values tuple.
    clusters : dict
        contains the cluster IDs as keys, and a centroid location, as well as a list of data points in the values tuple
    Returns
    -------
    clusters : dict
        contains the cluster IDs as keys, and the updated centroid location, as well as a list of data points in the values tuple

    """
    # loop over all clusters
    for key in clusters:
        # retrieve list of data points belonging to this cluster
        datapoints = clusters[key][1]
        
        # check that the cluster isn't empty
        if datapoints != []:
            
            # create an empty list that will contain the sum of all vectors of the data points
            vector_sum = [0]*len(data[datapoints[0]][2])
            
            # loop over all datapoints to calculate the sum of all vectors
            for x in datapoints:
                for i in range(len(vector_sum)):
                    vector_sum[i] += data[x][2][i]

            # compute the centroid vector 
            centroid = [x/len(datapoints) for x in vector_sum]
            
            # store the centroid location in the clusters dictionary 
            clusters[key][0] = centroid
    
    return clusters

def assign_datapoints(data: dict, clusters: dict, dist = 0): # 0 is euclidean, 1 is squared
    """
    Function that assigns each datapoint to the nearest cluster, based on the euclidean distance between the datapoint and each cluster centroid.
    
    Parameters
    ----------
    data : dict
        key-value pairs with data points as keys and values at the second place in the values tuple.
    clusters : dict
        contains the cluster IDs as keys, and their centroid location, as well as a list of data points in the values tuple
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    data : dict
        contains key-value pairs, with data points as keys and amongst the values are the updated cluster IDs
    clusters : dict
        contains the cluster IDs as keys, and as values the outdated centroid location, as well as an updated list of cell line indices

    """
    # empty list of data points
    for key in clusters:
        clusters[key][1] = []
    
    # loop over all data points in data
    for datapoint in data:
        centroids = []

        # loop over all clusters to calculate distance of data point to each centroid
        for cluster in clusters:
            if (dist == 0) or (dist == 1): distance = euclidean(data[datapoint][2], clusters[cluster][0], dist)
            
            # create list of cluster IDs with distances to it
            centroids.append([cluster, distance])

        # sort list of distances to centroids and choose the cluster ID and centroid belonging to the smallest distance
        nearest = sorted(centroids, key = lambda x:x[1])[0]
        
        # assign the cluster ID as attribute to the data point
        
        data[datapoint][3] = nearest[0]
        # assign the data point to the list of data points belonging to the centroid
        clusters[nearest[0]][1].append(datapoint)
        
    return data, clusters
        
def clustering(data: dict, k: int, max_iterations: int, dist = 0):
    """
    Function that performs the total clustering operation. It first generates a given amount of clusters.
    Then, for a given number of iterations, the cluster centroid positions are updated and datapoints are redistributed over the clusters.
    Meanwhile, it computes the Silhouette coefficient for each clustering, to check whether the clustering is still improving; if not the clustering is terminated.
    It also keeps count of the data types are in each cluster, and prints these stats for the final clustering.
    
    Parameters
    ----------
    data : dict
        contains key-value pairs, with data points as keys and amongst the values are the cluster IDs
    k : int
        number of clusters to be generated
    max_iterations : int
        maximal number of iterations, or times the clusters should be updated
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    data : dict
        contains key-value pairs, with data points as keys and amongst the values are the updated cluster IDs
    clusters : dict
        contains the cluster IDs as keys, and as values an updated list of data points, as well as the outdated centroid location
    silh : float
        Silhouette coefficient of the final clustering
        
    """
    # choose random centroids for the initial clustering
    clusters = initial_clustering(k, data);
    
    # initiate list of Silhouette coefficients
    S = [0]
    
    # iterate 
    for i in range(max_iterations):
        # assign datapoints to clusters
        data, clusters = assign_datapoints(data, clusters, dist = 0)
        
        # create new dictionary of non-zero clusters 
        non_zero_clusters = {}
        for key in clusters:
            if clusters[key][1] != []: non_zero_clusters[key] = clusters[key] 
        
        # create new clusters dictionairy with updated centroids, out of the non-zero clusters
        clusters = compute_centroid(data, non_zero_clusters)
        
        # calculate Silhouette coefficient and append to list
        S.append(silhouette_coeff(data, clusters, dist = 0))
    
        # stopping criterion
        if S[-1] == S[-2]: break
    
    # print final Silhouette score
    silh = S[-1]; #print('Silhouette coefficient:', silh, '\n')
    
    # keep count of data types
    for key in clusters:
        
        NB = 0
        BRCA = 0
        KIRC = 0
        COAD_READ = 0
        
        for datapoint in clusters[key][1]:
            if data[datapoint][1] == 'NB': NB += 1
            if data[datapoint][1] == 'BRCA': BRCA += 1
            if data[datapoint][1] == 'KIRC': KIRC += 1
            if data[datapoint][1] == 'COAD/READ': COAD_READ += 1
        
        # store stats in clusters dictionary:
        stats = {}
        stats['NB'] = NB; stats['BRCA'] = BRCA; stats['KIRC'] = KIRC; stats['COAD/READ'] = COAD_READ        
        clusters[key][2] = stats
        
    return data, clusters, silh