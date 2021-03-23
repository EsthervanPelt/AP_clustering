# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 09:51:15 2021

@author: 20183816
"""
def mean(x):
    return sum(x)/len(x)

def euclidean(u: list, v: list, squared = 1):
    if len(u) != len(v): print('Error: vectors not of the same length')
    
    elif len(u) == len(v):
        squared_euc = sum(map(lambda x, y: (x-y)**2, u, v))
        euc = squared_euc**(1/2)
        
        if squared == 0: return euc
        elif squared == 1: return squared_euc
    
def manhattan(u: list, v: list):
    if len(u) != len(v): print('Error: vectors not of the same length')
    
    elif len(u) == len(v):
        manh = sum(map(lambda x, y: abs(x-y), u, v))
        return manh

def correlation(u: list, v: list, dist = 1):
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
    NC = None
    for cluster in clusters.keys():
        if cluster != data[x][3]:
            MC = mean_cluster_weight(x, data, clusters[cluster], dist)
            
            if type(MC) != None:
                if (NC == None) or (NC > MC): NC = MC
    
    return NC

def silhouette_score(x: int, data: dict, clusters: dict, dist: 0): #0 = euclidean, 1 = squared euclidean
    cluster_x = data[x][3]
            
    W = mean_cluster_weight(x, data, clusters[cluster_x], dist)
    NC = neighbouring_cluster(x, data, clusters, dist) 

    if W < NC: s = 1-W/NC
    elif W == NC: s = 0
    elif W > NC: s = NC/W-1
    
    return s

def silhouette_coeff(data: dict, clusters: dict, dist = 0): #0 = euclidean, 1 = squared euclidean
    sum_s = sum(map(lambda x: silhouette_score(x, data, clusters, dist), data.keys()))
    S = sum_s / len(data.keys())
    
    return S

    
    