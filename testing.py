# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 09:22:44 2021

@author: 20183816
"""
import numpy as np
import time
import random
import math
import networkx as nx

def dot(u: list, v: list):
    time0 = time.time()
    #a = u@v
    time1 = time.time()
    b = u.T.dot(v)
    time2= time.time()
    c = sum([i*j for (i,j) in zip(u,v)])
    time3 = time.time()
    d =  sum(map(lambda x, y: x*y, u, v))
    time4 = time.time()
    print(#'@: '+time1-time0,a, 
          'numpy: '+str(time2-time1), b, '\n'+
          'zip: '+str(time3-time2), c, '\n'+
          'lambda: '+str(time4-time3), d)
    return d

def sqrt(x):
    time0 = time.time()
    a = x**(1/2)
    time1 = time.time()
    b = np.sqrt(x)
    time2 = time.time()
    c = math.sqrt(x)
    time3 = time.time()
    print('python: '+str(time1-time0), a, '\n'+
          'numpy: '+str(time2-time1), b, '\n'+ 
          'math: '+str(time3-time2), c)

u = np.random.rand(1000000,1)
v = np.random.rand(1000000,1)
d = dot(u,v)

sqrt(14234523454)


###### Test kmeans
from util import *
from kmeans import *
x = [0,1,2,5,5,5,8,9]
y = [5,3,5,1,8,11,1,0]
data = {}
for i in range(len(x)):
    data[i] = [None, None, [x[i],y[i]], None]
clusters = {}
clusters[1] = [[5,8],[]]
clusters[2] = [[5,11],[]]
clusters[3] = [[8,1],[]]

data, clusters = assign_datapoints(data, clusters, dist = 0) #1st iteration
clusters = compute_centroid(data, clusters)
data, clusters = assign_datapoints(data, clusters, dist = 0) #2nd iteration

###### Test working with graph of graphs
G = nx.MultiGraph()

for i in range(1,6):
    H = nx.Graph()
    H.add_nodes_from([i,i**2])

    G.add_node(H)

l = [0]+list(G.nodes)

G.add_edges_from([(l[1],l[2]), (l[1],l[3]), (l[1],l[5]), (l[2],l[4]), (l[2],l[3]), (l[3],l[4]), (l[4],l[5])])
list(G.nodes)[0].add_edge(1,0.1)
list(G.nodes)[1].add_edge(2,4)

for edge in G.edges(l[2]):
    if l[1] !=edge[1]:
        G.add_edge(l[1],edge[1])
        
print(list(G.nodes)[0].nodes, list(G.nodes)[0].edges)
print(list(G.nodes)[1].nodes, list(G.nodes)[1].edges)
l[1].add_nodes_from(l[2])
l[1].add_edges_from(l[2].edges)
print(list(G.nodes)[0].nodes, list(G.nodes)[0].edges)
G.remove_node(l[2])