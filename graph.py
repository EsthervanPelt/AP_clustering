# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:19:41 2021

@author: 20183816
"""
from util import half_correlation_matrix
from data import read_data
import matplotlib.pyplot as plt

# assign datafiles
fileMetadata = "GDSC_metadata.csv"
fileRMAExpression = "GDSC_RNA_expression.csv"

# read in data
data, gene_names = read_data(fileMetadata, fileRMAExpression);

corr_matrix = half_correlation_matrix(data,dist=0)

fc = []
y = 200
c = [x/(y-1) for x in range(y)]
for cc in c:
    f = 0
    n = len(corr_matrix)
    for row in corr_matrix:
        for i in row: 
            if i >= cc: f +=1
            
    f = f/(n*(n-1)/2)
    
    fc.append(f)

plt.figure()
plt.plot(c, fc)
plt.xlabel('c')
plt.ylabel('f(c)')
plt.title('fraction of nodes f(c) with correlation coefficient >= c')
plt.show()

cc = c[fc.index(list(filter(lambda i: i<=0.1, fc))[0])]

