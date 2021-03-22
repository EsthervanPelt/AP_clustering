# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 12:15:42 2021

@author: 20183816
"""
import random as rd
from util import *

def initial_clustering(k: int, data: dict):
    clusters = {}
    
    for cluster_ID in range(k):
        clusters[cluster_ID] = (None, []) #0 = centroid location, 1 = list of data points
    
    for key in data:
        cluster_ID = rd.randrange(k)
        key[3] = cluster_ID
        clusters[cluster_ID][1].append(key)
        
    return data, clusters

def compute_centroid(data: dict, clusters: dict):
    for key in clusters:
        datapoints = clusters[key][1]
        
        values = [0]*len(data(datapoints[0])[2])
        for x in datapoints:
            values = [sum(x) for x in zip(values, data[x][2])]
        centroid = [x/len(values) for x in values]
        
        clusters[key][0] = centroid
    
    return clusters

def assign_datapoints(data: dict, clusters: dict, dist = 0): # 0 is euclidean, 1 is squared
    for key in clusters:
        clusters[key][1] = []
    
    for datakey in data:
        centroids = []
        
        for clusterkey in clusters:
            if (dist == 0) or (dist == 1): distance = euclidean(data[datakey][2], clusters[clusterkey][0], dist)
            centroids.append((clusterkey, distance))
            
        nearest = centroids.min(key = lambda x:x[1])
        data[datakey][3] = nearest
        
        clusters[clusterkey][1].append(datakey)
        
    return data, clusters
        
        