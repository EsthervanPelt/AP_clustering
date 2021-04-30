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
        contains key-value pairs, with cell line indices as keys and cluster IDs as one of the values (this now contains Nones as placeholders)

    Returns
    -------
    clusters : dict
        contains the cluster IDs as keys, and a list of cell line indices as value

    """
    clusters = {}
    
    for cluster_ID in range(k):
        random_point = rd.choice(list(data.keys()))
        clusters[cluster_ID] = [data[random_point][2], []] #0 = centroid location, 1 = list of data points
        
    return clusters

def compute_centroid(data: dict, clusters: dict):
    """
    Function that computes the cluster centroid for each of the clusters. It uses the list of cell line indices contained in the clusters dictionairy,
    combined with the gene expression values of those cell lines that are contained within the data dictionairy. 

    Parameters
    ----------
    data : dict
        contains key-value pairs, with cell line indices as keys and gene expression values as one of the values
    clusters : dict
        contains the cluster IDs as keys, and as values a list of cell line indices, as well as the centroid location

    Returns
    -------
    clusters : dict
        contains the cluster IDs as keys, and as values a list of cell line indices, as well as the updated centroid location

    """
    for key in clusters:
        datapoints = clusters[key][1]
        
        if datapoints != []:
            values = [0]*len(data[datapoints[0]][2])
            for x in datapoints:
                for i in range(len(values)):
                    values[i] += data[x][2][i]

            centroid = [x/len(datapoints) for x in values]

            clusters[key][0] = centroid
    
    return clusters

def assign_datapoints(data: dict, clusters: dict, dist = 0): # 0 is euclidean, 1 is squared
    """
    Function that assigns each datapoint to the nearest cluster, based on the euclidean distance between the datapoint and each cluster centroid.
    
    Parameters
    ----------
    data : dict
        contains key-value pairs, with cell line indices as keys and gene expression values as one of the values
    clusters : dict
        contains the cluster IDs as keys, and as values a list of cell line indices, as well as the centroid location
    dist : boolean value
        indicates whether the normal or squared euclidean distance should be used.
        0 = normal, 1 = squared. The default is 0.

    Returns
    -------
    data : dict
        contains key-value pairs, with cell line indices as keys and amongst the values are the updated cluster IDs
    clusters : dict
        contains the cluster IDs as keys, and as values an updated list of cell line indices, as well as the outdated centroid location

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
    It also keeps count of how many of each cancer types are in each cluster, and prints these stats for the final clustering.
    
    Parameters
    ----------
    data : dict
        contains key-value pairs, with cell line indices as keys and gene expression values as one of the values
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
        contains key-value pairs, with cell line indices as keys and amongst the values are the updated cluster IDs
    clusters : dict
        contains the cluster IDs as keys, and as values an updated list of cell line indices, as well as the outdated centroid location
    silh : float
        Silhouette coefficient of the final clustering
        
    """
    # choose random centroids and assign datapoints
    clusters = initial_clustering(k, data);
    data, clusters = assign_datapoints(data, clusters, dist = 0)
    
    S = [0]
    # # iterate 
    for i in range(max_iterations):
        data, clusters = assign_datapoints(data, clusters, dist = 0)
        
        non_zero_clusters = {}
        for key in clusters:
            if clusters[key][1] != []: non_zero_clusters[key] = clusters[key] 
        
        clusters = compute_centroid(data, non_zero_clusters)
        
        S.append(silhouette_coeff(data, clusters, dist = 0))
    
        if S[-1] == S[-2]: break
    
    silh = S[-1]
    print('Silhouette coefficient:', silh, '\n')
    
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
            
        print('Cluster', key, '\n',
              'NB:       ', NB, '\n', 
              'BRCA:     ', BRCA, '\n',
              'KIRC:     ', KIRC, '\n',
              'COAD/READ:', COAD_READ)
    
    return data, clusters, silh