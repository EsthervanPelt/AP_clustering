# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 12:15:42 2021

@author: 20183816
"""
import random as rd
from util import euclidean
from util import silhouette_coeff

def initial_clustering(k: int, data: dict):
    clusters = {}
    
    for cluster_ID in range(k):
        random_point = rd.choice(list(data.keys()))
        clusters[cluster_ID] = [data[random_point][2], []] #0 = centroid location, 1 = list of data points
        
    return clusters

def compute_centroid(data: dict, clusters: dict):
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
    # choose random centroids and assign datapoints
    clusters = initial_clustering(k, data);
    data, clusters = assign_datapoints(data, clusters, dist = 0)
    
    S = [0]
    # # iterate 
    max_iterations = 20
    
    for i in range(max_iterations):
        data, clusters = assign_datapoints(data, clusters, dist = 0)
        
        non_zero_clusters = {}
        for key in clusters:
            if clusters[key][1] != []: non_zero_clusters[key] = clusters[key] 
        
        clusters = compute_centroid(data, non_zero_clusters)
        
        S.append(silhouette_coeff(data, clusters, dist = 0))
    
        if S[-1] == S[-2]: break
     
    print('Silhouette coefficient:', S[-1], '\n')
    
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
    
    return data, clusters, S[-1]