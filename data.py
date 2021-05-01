# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 10:14:03 2021

@author: 20183816
"""
import pandas as pd

def read_data(metadata: str, expression: str):
    """
    Function that reads in the gene expression values and metadata files (which are in csv format).
    The information is stored in a dictionary that contains the cell line indices, with the corresponding names, cancer types and gene expression values. 
    The cluster ID is not yet assigned, so this value contains Nones as a placeholder.

    Parameters
    ----------
    metadata : str
        filename of the csv-file containing metadata
    expression : str
        filename of the csv-file containing gene expression values

    Returns
    -------
    data : dict
        contains key-value pairs as specified above
    gene_names : list
        list of gene names, with indices corresponding to the set of gene expression values. 

    """
    # read in the data from the given files
    meta = pd.read_csv(metadata)
    genes = pd.read_csv(expression)
    
    # create empty data dictionary
    data = {}
    
    # fill the dictionary with information from the pandas matrices
    for index in meta.index:
        data[meta['COSMIC_ID'][index]] = [meta['name'][index],          # name of cel line
                                          meta['TCGA_label'][index],    # type of cancer
                                          genes.loc[genes['Unnamed: 0'] == meta['name'][index]].values[0,1:].tolist(), # gene expression values 
                                          None]                         # cluster ID
    
    # create a list containing the gene names
    gene_names = genes.columns.values[1:].tolist()
    
    return data, gene_names
 
# fileMetadata = "GDSC_metadata.csv"
# fileRMAExpression = "GDSC_RNA_expression.csv"

# data, gene_names = read_data(fileMetadata, fileRMAExpression)

