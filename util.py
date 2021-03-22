# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 09:51:15 2021

@author: 20183816
"""
def mean(x):
    return sum(x)/len(x)

def euclidean(u: list, v: list, squared = 1):
    squared_euc = sum(map(lambda x, y: (x*y)**2, u, v))
    euc = squared_euc**(1/2)
    
    if squared == 0: return euc
    elif squared == 1: return squared_euc
    
def manhattan(u: list, v: list):
    manh = sum(map(lambda x, y: abs(x-y), u, v))
    return manh

def correlation(u: list, v: list, dist = 1):
    mean_u = mean(u)
    mean_v = mean(v)
    
    nom = sum(map(lambda x,y: (x-mean_u)*(y-mean_v), u, v))
    den = sum(map(lambda x: (x-mean_u)**2))*sum(map(lambda y: (y-mean_v)**2))
    
    corr = nom/den 
    d_corr = 1 - abs(corr)
    
    if dist == 1: return d_corr
    elif dist == 0: return corr
    

    