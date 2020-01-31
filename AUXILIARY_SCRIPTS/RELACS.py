#!/usr/bin/env python

### functions for RELACS normalization ###

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.interpolate import interpn
import math
#import sys

def preprocess_deeptools(counts, log=False, merge_rep=False, rowSumFilt = 10, list_of_columns=[]):

    '''
    This function allows to sum counts between replicates or log-transform counts.
    To sum columns, provide a list of lists, where each sublist countain the indexes of the 
    columns that you want to merge.
    
    Ex. If you want to sum first and second columns together and third and fourth columns together, 
    the list_of_columns will be like: [[0,1],[2,3]]
    '''

    preprocessed_counts = {}
    
    chrom = list(counts[list(counts)[0]].values)
    start = list(counts[list(counts)[1]].values)
    end = list(counts[list(counts)[2]].values)
    
    counts.index = ['{}_{}_{}'.format(chrom[i],start[i],end[i]) for i in range(len(chrom))]
    counts = counts[list(counts)[3:]]
    # merge columns
    if merge_rep==True:
        if len(list_of_columns) > 0:
            for j in range(len(list_of_columns)):
                name = "+".join([list(counts)[i] for i in list_of_columns[j]])
                preprocessed_counts[name] = counts[[list(counts)[i] for i in list_of_columns[j]]].sum(axis=1).astype(np.int64)

            preprocessed_counts = pd.DataFrame(preprocessed_counts, index=counts.index)
        else:
            print("Please provide a list of columns that you want to merge.")
            return



    else:
        preprocessed_counts = counts[list(counts)].copy().astype(np.int64)
        
    preprocessed_counts = preprocessed_counts.loc[preprocessed_counts.sum(axis=1)>rowSumFilt]

    if log==True:
        
        preprocessed_counts = preprocessed_counts + 1

        preprocessed_counts = np.log2(preprocessed_counts)
        preprocessed_counts = preprocessed_counts.replace([np.inf,-np.inf],np.nan).dropna()


    return preprocessed_counts
    


def Normalize_T(counts, log=False, input_norm=True):
    
    '''
    This function takes as input a count dataframe and performs a CPM normalization 
    on the total library size for ChIP and Input. 
    Normalized counts can be output in log scale;
    if input_norm, function returns ChIP samples as ChIP/Input ratio.
    
    The dataframe is expected to have an even number of columns;
    - First half of columns are samples in the ChIP;
    - Second half of columns are the respective samples in the Input;
    - The order of samples in ChIP must be the same order as in Input.
    '''
    
    counts_norm = counts.copy()

    n_ = len(list(counts_norm))
    if n_%2 == 0:
        n_ = int(n_/2)
        tot_seq_Chip = counts_norm[list(counts_norm)[:n_]].sum().sum()
        tot_seq_I = counts_norm[list(counts_norm)[n_:]].sum().sum()
        
        #counts_norm[list(counts_norm)[:n_]] = (counts_norm[list(counts_norm)[:n_]].values / tot_seq_Chip) * 1000000
        #counts_norm[list(counts_norm)[n_:]] = (counts_norm[list(counts_norm)[n_:]].values / tot_seq_I) * 1000000

        if log == True:
            
            counts_norm = np.log2(counts_norm)
            counts_norm = counts_norm.replace([np.inf,-np.inf],np.nan).dropna()
            #print(counts_norm.shape)
            
            if input_norm == True:
                counts_norm[list(counts_norm)[:n_]] = counts_norm[list(counts_norm)[:n_]].values - counts_norm[list(counts_norm)[n_:]].values

        else:
            if input_norm == True:
                counts_norm[list(counts_norm)[:n_]] = np.log2(counts_norm[list(counts_norm)[:n_]].divide(counts_norm[list(counts_norm)[n_:]].values))
                counts_norm = counts_norm.replace([np.inf,-np.inf],np.nan).dropna()
    else:
        print("Your dataframe has an odd number of columns. Make sure that the first half of your DF\n\
        are ChIP, and the second half are the respective Input - The order of samples in the first half\n\
        must be the same as the order of samples in the second half!")
    
      
    return counts_norm




def Normalize_F(counts, cols, input_norm=True):
    
    '''
    This function takes as input a count dataframe. it first performs a CPM normalization 
    on the total library size for ChIP and Input. Then it computes sample specific scaling factors 
    on MT reads, which enforce sum of MT reads in ChIP to equal the sum of MT read in Input. 
    Normalized counts are output in log2 scale;
    if input_norm, function returns ChIP samples as ChIP/Input ratio.
    
    The dataframe is expected to have an even number of columns;
    - First half of columns are samples in the ChIP;
    - Second half of columns are the respective samples in the Input;
    - The order of samples in ChIP must be the same order as in Input.
    - the postional argument cols expects a list of lists, where the first element of each sublist is 
    a chip sample, and the second is its respective input.
    '''
    
    counts = Normalize_T(counts, log=False, input_norm=False)
    
    sf = []
    
    for i in range(len(cols)):
        
        chip_ = counts[cols[i][0]].loc[[i.split("_")[0].startswith("MT") for i in counts.index]].sum()
        inp_ = counts[cols[i][1]].loc[[i.split("_")[0].startswith("MT") for i in counts.index]].sum()
        sf_ = inp_/chip_
        sf.append(sf_)
                 
    counts[list(counts)[:len(sf)]] = counts[list(counts)[:len(sf)]] * sf

    counts = np.log2(counts)
    counts = counts.replace([np.inf,-np.inf],np.nan).dropna()
    
    if input_norm == True:
        counts[list(counts)[:len(sf)]] = counts[list(counts)[:len(sf)]].values - counts[list(counts)[len(sf):]].values
        
        

    return counts
    
  

def filter_counts(df_counts, rowsum = 0):
    
    '''
    Filters out genes whose rowsum is less or equal to a user defined threshold.
    By defalut, the function filters out genes that have zero reads aligned to them across all samples.
    '''
    
    df_filtered = df_counts.loc[df_counts.sum(axis=1) > rowsum]
    
    return df_filtered


def RLE(df_counts, rowSumFilt = 10):
    
    ''' 
    This function performs RLE normalization implemented in DESeq2 
    (this function outputs the same scaling factors as computed by DESeq2::estimateSizeFactorsForMatrix()
    and applies them to the input counts dataframe, returning a RLE normalized expression table).
    input: pandas DataFrame with gene ID as index (rows names) and sample names as columns (column names)
    output: 1. scaling factors
            2. RLE normalized count DataFrame
    '''
    
    # 1. Filter out untrascribed genes
    df_counts_filt = filter_counts(df_counts, rowsum = rowSumFilt)
    
    # 2. Take the ln of all values 
    df_log = np.log(df_counts_filt).replace([np.inf, -np.inf], np.nan).dropna()
    
    # 3. Compute the row-wise geometric average (geometric average of each gene expression across all samples)
    W_average = df_log.mean(axis = 1)
    
    # 4. Subtract the geometric mean from the log-transformed value of each gene
    df_sub = df_log.subtract(W_average, axis = "index")
    
    # 5. Calculate the median of the ratios for each sample
    median_ratios = df_sub.median(axis = 0)
    
    # 6. Convert the median from log space to normal number to get the scaling factors
    scaling_factors = math.e**(median_ratios)
    
    # 7. Divide the original read counts by the scaling factors
    df_norm = df_counts.divide(scaling_factors, axis = "columns")
    
    return scaling_factors, df_norm
    



def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins)
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter( x, y, c=z, **kwargs, s=1)

    return ax



def compare_samp(count, columns):

    n_ = len(columns)

    mt = count.loc[[i.split("_")[0].startswith("MT") for i in count.index]]

    fig,ax=plt.subplots(1,n_,dpi=300,figsize=(2+n_,2), sharey=True)

    for i in range(len(ax)):
        density_scatter(count[list(count)[n_:]].mean(axis=1), count[columns[i]], ax = ax[i], label="")
        ax[i].scatter(mt[list(count)[n_:]].mean(axis=1), mt[columns[i]], color="r",s=5, label="MT genes")
        ax[i].axhline(0,linestyle="--",linewidth=0.5,color='k')
        ax[i].set_title("{}".format(columns[i]),size=5)
    #ax[1].axhline(0,linestyle="--",linewidth=0.5,color='k')

    plt.tight_layout()

    return ax


def distpl_RELACS(counts, list_columns):

    # construct cmap
    my_cmap = sns.color_palette("Set2").as_hex()

    n_pl = len(list_columns)

    fig, ax = plt.subplots(n_pl,1,dpi=100,sharex=True,figsize=(5,8))

    for i in range(n_pl):
        idx=[list(counts)[i] for i in list_columns[i]]
        for k in range(len(idx)):
            sns.distplot(counts[idx[k]],ax=ax[i],color=my_cmap[k])
            sns.distplot(counts[idx[k]].loc[[i.split("_")[0].startswith("MT") for i in counts.index]], rug=True, hist=False, kde=False, ax=ax[i], rug_kws={"height":0.3},color=my_cmap[k])

    plt.tight_layout()

    return ax



def MA_samp(counts, list_columns, ax=None, reg=False):
    
    mt = counts.loc[[i.split("_")[0].startswith("MT") for i in counts.index]]
    
    if ax == None:
      pass
        # fig,ax = plt.subplots(dpi=400,figsize=(4,2.5))
    else:
        ax = ax
        
    lfc_2 =  counts[list_columns[1]] - counts[list_columns[0]]
    mean_cov = counts[list_columns].mean(axis=1)
    # density_scatter(mean_cov, lfc_2, ax=ax)
    #ax.scatter(counts[list_columns].loc[[i.split("_")[0].startswith("MT") for i in counts.index]].mean(axis=1),lfc_2.loc[[i.split("_")[0].startswith("MT") for i in counts.index]], s=3,color='r')
    if reg == True:
        sns.regplot(counts[list_columns].mean(axis=1),lfc_2,ax=ax, scatter=False, lowess=True, line_kws={'color':'r','linewidth':1,'linestyle':'--'})
    #ax.set_title("{}_vs_{}".format(list_columns[1],list_columns[0]), size=8)
    
    # ax.axhline(0,linestyle='--',linewidth=1, color='k')
    #plt.axhline(np.median(lfc_2),linestyle='--',linewidth=0.4, color='r')
    print(np.median(lfc_2))
    #plt.axvline(0,linestyle='--',linewidth=0.4, color='k')
    #ax.set_ylim((-8,8))
    plt.tight_layout()
    
    return ax, lfc_2, mean_cov
    
