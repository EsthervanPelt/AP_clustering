# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 10:14:03 2021

@author: 20183816
"""
import pandas as pd

def read_data(metadata: str, expression: str):
    meta = pd.read_csv(metadata)
    genes = pd.read_csv(expression)
    
    data = {}
    
    for index in meta.index:
        data[meta['COSMIC_ID'][index]] = (meta['name'][index],          # name of cel line
                                          meta['TCGA_label'][index],    # type of cancer
                                          genes.loc[values['Unnamed: 0'] == meta['name'][index]].values[0,1:].tolist(), # gene expression values 
                                          None)                         # centroid location
    
    gene_names = genes.columns.values[1:].tolist()
    return data, gene_names
 
fileMetadata = "GDSC_metadata.csv"
fileRMAExpression = "GDSC_RNA_expression.csv"

meta = pd.read_csv(fileMetadata)
values = pd.read_csv(fileRMAExpression)

data, gene_names = read_data(fileMetadata, fileRMAExpression)

