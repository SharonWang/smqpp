#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
import anndata
from scipy.stats import chi2
from scipy.sparse import issparse
import statsmodels.stats.multitest as multi
import statsmodels.api as sm
import warnings

try:
    from scanpy import logging as logg
except ImportError:
    pass

from .version import __version__

############# BELOW is data formating and readin #####################
def generate_feature_table(infile, outfile):
    '''
    Generate the feature table as input for read_in_files function
    
    Input
    -----
    infile: path to the gtf file that used for feature counting
    outfile: path to the output feature table in .tsv format
    
    Returns
    -----
    A feature table.tsv file as the input for read_in_files function
    
    '''
    
    savef = {}
    for line in open(infile):
        if '#' in line:
            continue
        val = line.split('\t')
        if val[2] == 'exon':
            ensemblID = re.search(r'.*gene_id "([^"]+).*', val[8]).group(1)
            if 'gene_name' not in val[8] or 'gene_biotype' not in val[8]:
                GN = ensemblID
                Type = 'Extra gene'
            else:
                GN = re.search(r'.*gene_name "([^"]+).*', val[8]).group(1)
                Type = re.search(r'.*gene_biotype "([^"]+).*', val[8]).group(1)
            savef[ensemblID] = [GN,Type]
    savef = pd.DataFrame.from_dict(savef, orient='index')
    savef.to_csv(outfile, index=True, header=False, sep='\t')

    
def reformat_meta(meta):
    '''
    Reformat metatable downloaded from google drive.
    This should be compatible with different versions of metadata spreadsheet.
    
    Input
    -----
    meta: metatable downloaded from google drive excel spreadsheet
    
    Returns
    -----
    Reformatted metatable
    
    '''
    
    if meta.shape[1] == 19:
        meta = meta.drop(meta.columns[[2,15]], axis=1)
    elif meta.shape[1] ==21:
        meta = meta.drop(meta.columns[[2,11,16,17]], axis=1)
    elif meta.shape[1] == 24:
        meta = meta.iloc[:,[0,1,2,3,4,5,6,9,7,8, 11, 10,21,22,17,16,15]]
    else:
        raise ValueError('To use this function, number of columns must be either 19, 21 or 24 depends versions on google drive.')
    #print(meta.columns)
    meta.columns = ['Gottgens_ID', 'Sequencing_identifier', 'Plate_number',
       'Position_in_96_well_plate_sorted', 'Position_in_96_well_plate_RNAseq',
       'FACs', 'Unique_sample_descriptor', 'Details', 'Cell_type_general',
       'Cell_type_subtype', 'Owner', 'Species', 'Sequencing_index',
       'Sequencing_facility_index', 'Average_pooled_library_length_bp', 'Pool_size',
       'Number_of_lanes']
    #print(meta.columns)
    return(meta)

    
def diff_list(list1, list2):
    c = set(list1).union(set(list2))  
    d = set(list1).intersection(set(list2)) 
    return list(c - d)


def read_in_files(Indir, ftable ,method = 'FeatureCount'):
    '''
    Read in STAR aligned files. 
    
    Input
    -----
    Indir: FeatureCount output folder path or HTSeqcount output file path
    ftable: Gene feature table, ftable must have two columns with ['Gene Name', 'Gene Type'] with Ensembl ID as index
    method: Counting method, can be either 'FeatureCount' or 'HTSeqcount', default: FeatureCount
    
    Returns
    -----
    An anndata object
    
    '''
    
    if method == 'FeatureCount':
        # first tidy up X
        X = pd.read_csv(Indir+'/fcounts/fcounts.txt', delimiter='\t', index_col=0, skiprows=1)
        GT = X.iloc[['ERCC-' not in x for x in X.index],0:5].copy()
        X = X.iloc[:,5:].transpose().copy()
        X['obs_names'] = [re.search('(SLX-\d+\.\w\d+_\w\d+)', x).group(0) for x in X.index]
        X = X.groupby('obs_names').sum()
        ERCC = X.iloc[:,['ERCC-' in x for x in X.columns]].copy()
        X = X.iloc[:,['ERCC-' not in x for x in X.columns]].copy()

        # Then work on the QC
        QC = pd.read_csv(Indir+'/fcounts/fcounts.txt.summary', delimiter='\t', index_col=0).transpose()
        QC['obs_names'] = [re.search('(SLX-\d+\.\w\d+_\w\d+)', x).group(0) for x in QC.index]
        QC = QC.groupby('obs_names').sum()
        QC = QC.iloc[:,(np.sum(QC, axis=0)!=0).values].copy()
        QC.columns = ['QC_'+x for x in QC.columns]
        QC = QC[[x for x in QC.columns if 'Unassigned' in x]]
    elif method == 'HTSeqcount':
        X = pd.read_csv(Indir, delimiter='\t', index_col=0).transpose()
        X['obs_names'] = X.index
        X = X.groupby('obs_names').sum()
        ERCC = X.iloc[:,['ERCC-' in x for x in X.columns]].copy()
        QC = X.iloc[:,['__' in x for x in X.columns]].copy()
        QC.columns = [x.replace('__','') for x in QC.columns]
        QC.columns = ['QC_'+x for x in QC.columns]
        X = X.iloc[:,['ERCC-' not in x and '__' not in x for x in X.columns]].copy()
    else:
        raise ValueError('method can only be either FeatureCount or HTSeqcount.')
    
    # Combine feature table
    # ftable = pd.read_csv(ftable_loc, delimiter='\t', index_col=0, header=None)
    ftable = ftable.iloc[['ERCC-' not in x for x in ftable.index],:].copy()
    difG = diff_list(list(X.columns), list(ftable.index))
    if difG:
        raise ValueError(f'Gene names in Feature table do not match the ones in count table. Different genes: {difG}')
    ftable.columns = ['Gene Name', 'Gene Type']
    ftable = ftable.loc[X.columns,:].copy()
    if method == 'FeatureCount':
        FeatureTable = pd.concat([ftable,GT], axis=1)
    else:
        FeatureTable = ftable
    print('Count table shape: '+str(X.shape))
    print('Feature table shape:' + str(FeatureTable.shape))
    
    # New write in anndata frame
    adata = anndata.AnnData(X=X, var=FeatureTable, obs=QC)
    
    # Here change ERCC to array as obsm only accepts arrays
    ERCC = np.array(ERCC)
    if ERCC.size == 0:
        warnings.warn('ERCC empty!')
    else:
        adata.obsm['ERCC'] = ERCC
    adata.var['Ensembl_ID'] = adata.var_names
    adata.var_names = adata.var['Gene Name']
    adata.var_names_make_unique()
    return adata


############## BELOW is QC #########################
cutoff = {}
cutoff['nMapped (log10)'] = np.log10(2*(10**5))
cutoff['nNuclear (log10)'] = 0
cutoff['fGenes:nTotal'] = 0.2
cutoff['nHCGenes'] = 4000
cutoff['mito:nGenes'] = 0.2
cutoff['nERCC:nMapped'] = 0.2

def smartseq_qc(adata, cutoff=cutoff,
            MTpattern = 'mt-', ncols=4, figsize=(10,7), s=10, title=None, save=None):
    '''
    Do a bglab equivalent quality control.
    
    Input
    -----
    adata: anndata object
    cutoff: QC cutoffs (dictionary), 6 in total, default: 
        {'nMapped (log10)': np.log10(2*(10**5)),
        'nNuclear (log10)': 0,
        'fGenes:nTotal': 0.2,
        'nHCGenes': 4000,
        'mito:nGenes': 0.2,
        'nERCC:nMapped': 0.2
        }
    MTpattern: mitochrondria gene pattern, default 'mt-' for mm10
    ncols: number of columns for plotting, default: 4
    figsize: size of the figure, default: (10,7)
    s: point size, default: 10
    title: tilte of the figure, default: None
    save: if save is not None, figure will be saved as named in save, default: None
    
    Returns
    -----
    An anndata object with cells passed QC
    n_counts, n_genes and percent_mito added in .obs
    
    '''
    ## don't know why make_unique did not work sometimes, so do this again
    adata.var_names_make_unique()
    mito_genes = [name for name in adata.var_names if name.startswith(MTpattern)]
    if not mito_genes:
        raise ValueError('Pls enter correct MTpattern')
    else:
        print('mito_genes: '+str(mito_genes))
    mitoCNT = np.sum(adata[:,mito_genes].X, axis=1).copy()
    nuclearCNT = np.sum(adata[:,~np.in1d(adata.var_names,mito_genes)].X, axis=1).copy()
    if 'ERCC' not in adata.obsm_keys():
        erccCNT = np.zeros(adata.shape[0])
    else:
        erccCNT = np.sum(adata.obsm['ERCC'], axis=1)
    qcNames = [x for x in adata.obs_keys() if 'QC' in x]
    if not qcNames:
        qcCNT = np.zeros(adata.shape[0])
    else:
        qcCNT = np.sum(adata.obs[qcNames], axis=1).values
    nTotal = mitoCNT + nuclearCNT + erccCNT + qcCNT
    nMapped = mitoCNT + nuclearCNT + erccCNT
    nGenes = mitoCNT + nuclearCNT
    nHCGenes = np.sum(adata[:,~np.in1d(adata.var_names,mito_genes)].X.T*(10**6)/(nuclearCNT+1) > 10, axis=0)

    QCdata = {}
    QCdata['nMapped (log10)'] = np.log10(nMapped+1)
    QCdata['nNuclear (log10)'] = np.log10(nuclearCNT+1)
    QCdata['fGenes:nTotal'] = nGenes/(nTotal+1)
    QCdata['nHCGenes'] = nHCGenes
    QCdata['mito:nGenes'] = mitoCNT/(nGenes+1)
    QCdata['nERCC:nMapped'] = erccCNT/(nMapped+1)
    QCdata['nNuclear:nMapped'] = nuclearCNT/(nMapped+1)
    for qcIndex in qcNames:
        QCdata[qcIndex.replace('_Unassigned', '')+':nTotal'] = adata.obs[qcIndex]/(nTotal+1)
    
    # cells failed QC
    compara = {}
    compara['nMapped (log10)'] = '<'
    compara['nNuclear (log10)'] = '<'
    compara['fGenes:nTotal'] = '<'
    compara['nHCGenes'] = '<'
    compara['mito:nGenes'] = '>'
    compara['nERCC:nMapped'] = '>'
    failed = []
    for k in cutoff.keys():
        if compara[k] == '<':
            failed.append(QCdata[k] < cutoff[k])
        else:
            failed.append(QCdata[k] > cutoff[k])
    failed = np.vstack(failed)
    failed_idx = np.sum(failed, axis=0)>0
    print('Number of passed cells: '+str(sum(np.sum(failed, axis=0)==0)))
    print('Number of failed cells: '+str(sum(failed_idx)))
    
    # plotting
    nrows = int(np.ceil(len(QCdata.keys())/ncols))
    #print(nrows)
    fig, ax = plt.subplots(nrows,ncols, figsize=figsize)
    for i in range(len(QCdata.keys())):
        colidx = i%ncols
        rowidx = np.floor(i/ncols).astype(int)
        ax[rowidx, colidx].scatter(nTotal, QCdata[list(QCdata.keys())[i]], s=s, color='black')
        ax[rowidx, colidx].scatter(nTotal[failed_idx], QCdata[list(QCdata.keys())[i]][failed_idx], s=s, color='red')
        if list(QCdata.keys())[i] in cutoff.keys():
            ax[rowidx, colidx].axhline(y=cutoff[list(QCdata.keys())[i]], color='orange', linestyle='dashed')
        #ax[rowidx, colidx].set_yscale('log',basey=10)
        ax[rowidx, colidx].set_ylabel(list(QCdata.keys())[i])
        ax[rowidx, colidx].grid(False)
    fig.text(0.5, -0.03, 'nTotal', ha='center')
    fig.suptitle(title)
    plt.tight_layout(pad=1)
    fig.subplots_adjust(top=0.88)
    
    if save is not None:
        plt.savefig(save)
        
    adata.obs['n_counts'] = nGenes
    adata.obs['percent_mito'] = mitoCNT/(nGenes+1)
    adata.obs['n_genes'] = np.sum(adata.X > 0, axis=1)
    
    return adata[~failed_idx,:].copy()


############## BELOW is DESeq2 normalisation ######
def fexp_genes(x):
    '''
    Get genes expressed in all cells
    
    '''
    return np.sum(x==0, axis=0)==0
    
def exp_genes(x):
    '''
    Get genes expressed at least in 1 cell
    
    '''
    return np.sum(x>0, axis=0)>0

def est_size_factor(x, method='ExpAllC'):
    '''
    Size factor estimation using DESeq2 normalisation method
    
    Input
    -----
    x: raw count matrix, rows: cells, columns: genes
    method: ['ExpAllC', 'Exp'], ExpAllC - Get genes expressed in all cells, Exp - Get genes expressed at least in 1 cell. Default: ExpAllC
    
    Returns
    -----
    Estimated size factor
    
    '''
    
    if method == 'ExpAllC':
        x = x[:,fexp_genes(x)].copy()
    elif method == 'Exp':
        x = x[:,exp_genes(x)].copy()
    else:
        raise Exception('Method needs to be either ExpAllC or Exp')
        
    print('Filtered matrix shape: '+ str(x.shape))
    loggeomeans = np.mean(np.log(x), axis=0)
    print('Number of valid means:' + str(sum(np.isfinite(loggeomeans))))
    sf = np.exp(np.median((np.log(x)-loggeomeans)[:,np.isfinite(loggeomeans)], axis=1))
    return sf

def normalise_data(adata, reCalSF=True, method='ExpAllC', copy=False):
    '''
    Normalisation using size factors estimated by DESeq2 method
    
    Input
    -----
    adata: an anndata object
    method: ['ExpAllC', 'Exp'], ExpAllC - Get genes expressed in all cells, Exp - Get genes expressed at least in 1 cell. Default: ExpAllC
    reCalSF: If recalculate the size factors, default: True, If False, 'sf_genes' and 'sf_ercc' in .obs will be used.
    copy: if copy anndata into new object, default: False
    
    Returns
    -----
    Updated anndata object
    sf_gene and sf_ercc are added in .obs
    ERCC_norm is added in .obsm
    
    '''
    
    if reCalSF:
        print('Calculate SF for genes:')
        sf_genes = est_size_factor(adata.X, method=method)
        adata.obs['sf_gene'] = sf_genes
        if 'ERCC' in adata.obsm_keys():
            if adata.obsm['ERCC'].size !=0:
                print('Calculate SF for erccs:')
                sf_ercc = est_size_factor(adata.obsm['ERCC'], method=method)
                adata.obs['sf_ercc'] = sf_ercc
    else:
        if 'sf_gene' not in adata.obs_keys():
            raise ValueError('sf_gene is not found in .obs, please set reCalSF=True.')
    
    adata.X = np.log1p(adata.X/adata.obs['sf_gene'].values[:,None])
    if 'ERCC' in adata.obsm_keys():
        if adata.obsm['ERCC'].size !=0:
            adata.obsm['ERCC_norm'] = np.log1p(adata.obsm['ERCC']/adata.obs['sf_ercc'].values[:,None])
    if copy:
        return adata.copy()

############## BELOW is Quantile normalisation ########
from scipy.stats import rankdata
def quantile_norm(X):
    '''
    Quantile normalisation
    
    Input
    -----
    X: adata.X, [cells, genes]. Note if adata.X is a sparse matrix, please first do X = adata.X.toarray()

    Returns
    -----
    Normalised counts, can be inserted back to adata.X
    
    '''
    
    quantiles = np.mean(np.sort(X, axis=1), axis=0)
    ranks_min = np.apply_along_axis(rankdata, 1, X, 'min')
    rank_min_indices = ranks_min.astype(int)-1
    ranks_max = np.apply_along_axis(rankdata, 1, X, 'max')
    rank_max_indices = ranks_max.astype(int)-1
    Xn_min = quantiles[rank_min_indices]
    Xn_max = quantiles[rank_max_indices]
    Xn = (Xn_min+Xn_max)/2
    Xn[X==0] = 0
    return(Xn)

def quantile_norm_log(X):
    '''
    Log quantile normalisation
    
    Input
    -----
    X: adata.X, [cells, genes]. Note if adata.X is a sparse matrix, please first do X = adata.X.toarray()

    Returns
    -----
    Normalised counts in log scale, can be inserted back to adata.X
    
    '''
        
    logX = np.log1p(X)
    logXn = quantile_norm(logX)
    return(logXn)

############## BELOW is downsampling normalisation #####
def downsampling(X, min_lib_size, seed=0):
    '''
    Downsampling normlisation for each cell, randomly generate a number using binomial distribution,
    with probability equal to the specific capture efficiency.
    
    Input
    -----
    X: raw UMI counts for each cell. Note if adata.X is a sparse matrix, please first do X = adata.X.toarray()
    min_lib_size: minimum library size among all cells
    seed: random seed, default: 0
    
    Returns
    -----
    Normalised counts for each cell
    
    '''
        
    np.random.seed(seed=seed)
    prob = min_lib_size/np.sum(X)
    return(np.array([np.random.binomial(x, prob, 1) for x in X]).flatten())

def downsampling_norm(X, seed=0):
    '''
    Downsampling normlisation for each cell, randomly generate a number using binomial distribution,
    with probability equal to the specific capture efficiency.
    
    Input
    -----
    X: adata.X, [cells, genes]. Note if adata.X is a sparse matrix, please first do X = adata.X.toarray()
    seed: random seed, default: 0
    
    Returns
    -----
    Normalised counts, can be inserted back to adata.X
    
    '''
        
    min_lib_size = np.min(np.sum(X, axis=1))
    return(np.apply_along_axis(downsampling, 1, X, min_lib_size, seed))

############## Below is HVG selection ##################
def tech_var(adata, useERCC=True, cvThresh=.3, quant=.8, minBiolDisp=.5**2, 
           fdr=.1, meanForFit=None, copy=False):
    '''
    Calculate highly variable genes (HVG) using Brennecke et. al method
    
    Input
    -----
    adata: an anndata object
    useERCC: if to use ERCC for HVG calculation, default: True
    cvThresh: The CV threshold to use for choosing cells to fit. See quant. Default: 0.3
    quant: The quantile of data to fit so that linear model is not fitted to plateau at lower gene expression values. Default: 0.8
    minBiolDisp: Assumed biological dispersion. Default: 0.25
    fdr: False discovery rate for chi-squared test. Default: 0.1
    meanForFit: Provide a minimum mean for fit. Default: None
    copy: if copy into a new anndata object. Default: False
    
    Returns
    -----
    Updated anndata object
    varGenes is added in .uns
    .uns['varGenes']['parameters']: parameters for modal fitting
    .uns['varGenes']['genes']: mean, cv2 and HVG index for genes
    .uns['varGenes']['ercc']: mean, cv2 for ERCCs
    '''
    
    if 'sf_gene' not in adata.obs_keys():
        print("sf_gene is not found, redoing normalisation")
        normalise_data(adata)
    
    data = np.exp(adata.X)-1
    aMean = np.mean(data, axis=0)
    aStd = np.std(data, axis=0)
    cv2a = (aStd/aMean)**2
    if useERCC:
        if 'ERCC' not in adata.obsm_keys():
            raise ValueError('ERCC does not exist, check data.')
        ercc_data = np.exp(adata.obsm['ERCC_norm'])-1
        sMean = np.mean(ercc_data, axis=0)
        sStd = np.std(ercc_data, axis=0)
        cv2s = (sStd/sMean)**2
    else:
        sMean = aMean
        cv2s = cv2a
    if meanForFit is None:
        meanForFit = np.quantile(sMean[cv2s>cvThresh], quant)
    print('MeanForFit: ', str(meanForFit))
    useForFit = (sMean>=meanForFit)
    print(np.sum(useForFit))
    
    a1tilde = 1/sMean[useForFit]
    x = sm.add_constant(a1tilde, prepend=False)
    y = cv2s[useForFit]
    link_func = sm.genmod.families.links.Identity()
    fit = sm.GLM(y, x, family=sm.families.Gamma(link=link_func)).fit()

    a0 = fit.params[1]
    a1t = fit.params[0]
    df = data.shape[0]-1
    m = data.shape[0]
    xi = None
    
    if useERCC:
        xi = np.mean(1/adata.obs['sf_ercc'])
        psi = xi + (a1t - xi)*np.mean(adata.obs['sf_ercc']/adata.obs['sf_gene'])
        cv2th = a0 + minBiolDisp + a0*minBiolDisp
        testDenom = (aMean*psi + cv2th*aMean**2)/(1+cv2th/m)
        pA = 1 - chi2.cdf((aStd**2)*df/testDenom, df=df)
    else:
        psi = a1t
        chi2_values = df * cv2s / (psi / sMean + a0)
        pA = 1 - chi2.cdf(chi2_values ,df=df)
    
    pA[np.isnan(pA)] = 1
    _, padj, _, _ = multi.multipletests(pA, method='fdr_bh')
    
    highVarGenes = adata.var_names[padj < fdr]
    print('Length of HVGs: '+ str(len(highVarGenes)))
    
    adata.uns['varGenes'] = {}
    adata.uns['varGenes']['parameters'] = {}
    adata.uns['varGenes']['genes'] = {}
    adata.uns['varGenes']['ercc'] = {}
    
    adata.uns['varGenes']['parameters']['minBiolDisp'] = minBiolDisp
    adata.uns['varGenes']['parameters']['a1tilde'] = a1t
    adata.uns['varGenes']['parameters']['a0'] = a0
    adata.uns['varGenes']['parameters']['psi'] = psi
    adata.uns['varGenes']['parameters']['xi'] = xi
    adata.uns['varGenes']['parameters']['df'] = df
    adata.uns['varGenes']['parameters']['useERCC'] = useERCC
    adata.uns['varGenes']['parameters']['meanForFit'] = meanForFit
    adata.uns['varGenes']['parameters']['useForFit'] = np.sum(useForFit)
    adata.uns['varGenes']['parameters']['cvThresh'] = cvThresh
    adata.uns['varGenes']['parameters']['quant'] = quant
    
    adata.uns['varGenes']['genes']['mean'] = aMean
    adata.uns['varGenes']['genes']['cv2'] = cv2a
    adata.uns['varGenes']['genes']['highVar'] = (padj < fdr)
    adata.uns['varGenes']['ercc']['mean'] = sMean
    adata.uns['varGenes']['ercc']['cv2'] = cv2s
    
    if copy:
        return adata.copy()

def plot_tech_var(adata, s=10, save=None):
    '''
    Plot the results for highly variable gene fitting.
    
    Input
    -----
    adata: an anndata object
    s: point size, default: 10
    save: if save is not None, figure will be saved as named in save, default: None
    
    Returns
    -----
    A figure showing HVG selection
    
    '''
    
    minBiolDisp = adata.uns['varGenes']['parameters']['minBiolDisp']
    a1t = adata.uns['varGenes']['parameters']['a1tilde']
    a0 = adata.uns['varGenes']['parameters']['a0']
    psi = adata.uns['varGenes']['parameters']['psi']
    xi = adata.uns['varGenes']['parameters']['xi'] 
    df = adata.uns['varGenes']['parameters']['df']
    useERCC = adata.uns['varGenes']['parameters']['useERCC'] 
    
    g_df = pd.DataFrame(adata.uns['varGenes']['genes'])
    g_df.index = adata.var_names
    g_df = g_df.loc[g_df['mean']>0,]
    
    e_df = pd.DataFrame(adata.uns['varGenes']['ercc'])
    e_df = e_df.loc[e_df['mean']>0,]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.scatter(g_df.loc[~g_df['highVar'],:]['mean'], g_df.loc[~g_df['highVar'],:]['cv2'], color='grey', s=s)
    ax.set_yscale('log',base=10)
    ax.set_xscale('log',base=10)
    ax.set_ylabel(r'$\sigma^{2}/\mu^{2}$')
    ax.set_xlabel(r'$\mu$')
    ax.grid(False)

    plt.scatter(g_df.loc[g_df['highVar'],:]['mean'], g_df.loc[g_df['highVar'],:]['cv2'], color='red', s=s)
    
    xlim = np.log10(ax.get_xlim())
    ylim = ax.get_ylim()
    xpred = 10**(np.linspace(start=xlim[0], stop=xlim[1], num=1000))
    if useERCC:
        plt.scatter(e_df['mean'], e_df['cv2'], s=s, color='blue')
        y1 = a1t/xpred + a0
        y2 = psi/xpred + a0 + minBiolDisp
        plt.plot(xpred[y1<ylim[1]], y1[y1<ylim[1]], color='red')
        plt.plot(xpred[y2<ylim[1]], y2[y2<ylim[1]], color='red', linestyle='dashed')
    else:
        y1 = psi/xpred + a0
        plt.plot(xpred[y1<ylim[1]], y1[y1<ylim[1]], color='red', linestyle='dashed')
    
    if save is not None:
        plt.savefig(save)

def detect_outlier_cells(adata, aMeanQ = 0.95, cv2aQ = 0.8, outQ = 0.8, s=10):
    '''
    Outlier cell detection
    
    Input
    -----
    adata: an anndata object
    aMeanQ: Quantile of mu value (x axis) to select. Default: 0.95
    cv2aQ: Quantile of cv2a value (y axis) to select. Default: 0.8
    outQ: percentage of selected genes to be consider as outlier. Default: 0.8
    s: point size. Default: 10
    
    Returns
    -----
    MA plot
    up-regulated gene list, down-regulated gene list, full DE table
    
    '''
    if 'varGenes' not in adata.uns_keys():
        raise ValueError('Please first do highly variable genes selection using tech_var function')
    
    aMean = adata.uns['varGenes']['genes']['mean']
    cv2a = adata.uns['varGenes']['genes']['cv2']
    idx = (aMean > np.quantile(aMean, aMeanQ)) & (cv2a > np.quantile(cv2a, cv2aQ))
    GL = adata.var_names[idx]
    print('Number of selected Genes: ' + str(len(GL)))
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.scatter(aMean, cv2a, color='grey', s=s)
    ax.set_yscale('log',basey=10)
    ax.set_xscale('log',basex=10)
    ax.set_ylabel(r'$\sigma^{2}/\mu^{2}$')
    ax.set_xlabel(r'$\mu$')
    ax.grid(False)
    plt.scatter(aMean[idx], cv2a[idx], color='red', s=s)
    
    GLData = adata[:,GL].X.copy()
    gmean = np.mean(GLData, axis = 0)
    gstd = np.std(GLData, axis = 0)
    outlierC = adata.obs_names[np.sum(GLData > gmean+gstd, axis=1) > len(GL)*0.8].values
    print('Number of outlier cell: '+ str(len(outlierC)))
    print('Outlier cells: '+ str(outlierC))
    return outlierC
        
        
########## BELOW is for DE analysis #####################
def plot_ma(adata, unsName='rank_genes_groups', cidx=0, Cells = None, save=False, padj_cutoff=0.05, logFC_cutoff=1, 
            exp_cutoff=-6, s=1):
    '''
    MA plot for differential expression analysis and select the significant genes with high confidence
    
    Input
    -----
    adata: an anndata object
    unsName: key_added from rank_genes_groups. Default: rank_genes_groups
    cidx: index of columns from rank_genes_groups results. Default: 0
    Cells: subset cells for the plot. Default: None
    save: if save is not None, figure will be saved as named in save. Default: None
    padj_cutoff: cutoff for adjusted p value. Default: 0.05
    logFC_cutoff: cutoff for log fold change. Default: 1
    exp_cutoff: cutoff for log mean expression. Default: -6
    s: point size. Default: 1
    
    Returns
    -----
    MA plot
    up-regulated gene list, down-regulated gene list, full DE table
    
    '''
    
    if Cells is not None:
        adata_sub = adata[Cells,:]
    else:
        adata_sub = adata
    print(adata_sub.shape)
    gnames = pd.DataFrame(adata.uns[unsName]['names']).iloc[:,cidx]
    logFC = pd.DataFrame(adata.uns[unsName]['logfoldchanges']).iloc[:,cidx]
    pvals = pd.DataFrame(adata.uns[unsName]['pvals']).iloc[:,cidx]
    padj = pd.DataFrame(adata.uns[unsName]['pvals_adj']).iloc[:,cidx]
    adata_sub = adata_sub.raw[:, gnames].X
    print(adata_sub.shape)
    if issparse(adata_sub):
        adata_sub = adata_sub.toarray()
    normExp = np.mean(np.exp(adata_sub)-1, axis=0)
    del adata_sub

    abs_logFC = logFC.copy()
    abs_logFC[abs_logFC > 4] = 4
    abs_logFC[abs_logFC < -4] = -4

    logExp = np.log2(normExp)
    idx = (padj < padj_cutoff) & (abs(abs_logFC) > logFC_cutoff)
    upidx = (padj < padj_cutoff) & (abs_logFC > logFC_cutoff) & (logExp > exp_cutoff)
    downidx = (padj < padj_cutoff) & (abs_logFC < -logFC_cutoff) & (logExp > exp_cutoff)
    print('upRegulated gene: '+str(sum(upidx)))
    print('downRegulated gene: '+str(sum(downidx)))
    
    fig = plt.figure()
    plt.scatter(x=logExp, y=abs_logFC, s=s)
    plt.scatter(x=logExp[idx & (logExp > exp_cutoff)], y=abs_logFC[idx & (logExp > exp_cutoff)], c='red',s=s)
    plt.axhline(y=0, color='black')
    plt.axhline(y=logFC_cutoff, color='grey', linestyle = '--')
    plt.axhline(y=-logFC_cutoff, color='grey', linestyle = '--')
    plt.axvline(x=exp_cutoff, color='grey', linestyle = '--')
    plt.xlabel('log2 Mean Exp')
    plt.ylabel('log2 Fold Change')
    plt.grid(b=None)
    plt.show()
    if save:
        fig.savefig(save)
    
    Ftable = pd.DataFrame(np.column_stack([gnames, logExp, logFC, pvals, padj]), columns=['GN','logMeanExp', 'logFC', 'pvals', 'padj'])
    return gnames[upidx], gnames[downidx], Ftable


###### BELOW is for pseudotime analysis ##################
import patsy
import numpy as np
from scipy.stats import chi2
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
from scipy.signal import gaussian

def ns(pt, df=3):
    '''
    Fitting the gene expression to a smooth, nonlinear function of pseudotime with natural spline 
    
    Input
    -----
    pt: pseudotime
    df: degree of freedom, default: 3
    
    Returns
    -----
    Natural spline fit
    
    '''
    
    return patsy.dmatrix('cr(x, df=3) -  1', {'x': pt})

def likelihood_ratio_test(X_alt, y, X_null=None):
    '''
    Computer likelihood ratio test between two models: full model and the reduced model
    if X_null is None: the reduced model will be only trained on the intercept. Note that
    reduced model must be a subset of full model -- it can not contain features that are 
    not in full model.
    
    Input
    -----
    X_alt: full model design matrix
    y: expression values for each gene
    X_null: reduced model design matrix
    
    Returns
    -----
    p-value, which can be used to accept or reject the null hypothesis
    
    '''
    
    y = np.array(y)
    X_alt = np.array(X_alt)
    X_alt = sm.add_constant(X_alt)
    
    if X_null is not None:
        X_null = np.array(X_null)
        X_null = sm.add_constant(X_null)
        
        if X_null.shape[1] >= X_alt.shape[1]:
            raise ValueError("Alternate features must have more features than null features")
        
        df = X_alt.shape[1] - X_null.shape[1]
    else:
        X_null = np.repeat(1,X_alt.shape[0])
        df = X_alt.shape[1]
    
    fit_alt = sm.OLS(y, X_alt).fit()
    fit_null = sm.OLS(y, X_null).fit()
    
    llf_alt = fit_alt.llf
    llf_null = fit_null.llf

    G = 2 * (llf_alt - llf_null)
    p_value = chi2.sf(G, df)

    return p_value

def GeneExp_LLR_test(adata, alt_obs, useHVG=True, null_obs=None, ns_df=3):
    '''
    Perform likelihood ratio test to detect genes that are differential expression along
    pseudotime. In order to do this, dpt_pseudotime must be in the .obs. 
    
    This can also be applied for general model comparison if dpt_pseudotime does not exist.
    
    Please NOTE that, the input .raw must be log and normalised expression values as the assumption
    is that the log norm exp should follow guassian distribution. If the raw counts are applied, this
    model cannot be used and more complicated generalised linear models should be used for fitting.
    
    Input
    -----
    adata: adata object
    alt_obs: .obs terms that will be considered in the full model, dpt_pseudotime must be in .obs to test for DE along PT 
    useHVG: only use highly variable gene (HVG), default: True
    null_obs: .obs terms that will be considered in the reduced model
    ns_df: degree of freedom for smoothing the PT, default: 3
    
    Returns
    -----
    Gene table with p values (pval) and adjusted p values (padj)
    
    '''
    
    if useHVG:
        yall = adata.raw[:,adata.var_names].X
        GN = adata.var_names
    else:
        yall = adata.raw.X
        GN = adata.raw.var_names
    Ngenes = yall.shape[1]
    Ncells = yall.shape[0]
    alt_obs = np.array(alt_obs)
    null_obs = np.array(null_obs)
    if 'dpt_pseudotime' in alt_obs:
        pt_design = ns(adata.obs['dpt_pseudotime'], df=ns_df)
        if len(alt_obs) > 1:
            alt_obs = alt_obs[alt_obs != 'dpt_pseudotime']
            Naltobs = len(alt_obs)
            alt_obs = adata.obs[alt_obs].to_numpy()
            alt_obs = alt_obs.reshape(Ncells, Naltobs)
            alt_obs = np.concatenate((alt_obs, pt_design), axis=1)
        else:
            alt_obs = pt_design
    else:
        alt_obs = adata.obs[alt_obs].to_numpy()
    alt_obs = alt_obs.astype(float)
    
    if null_obs is not None:
        if 'dpt_pseudotime' in null_obs:
            pt_design = ns(adata.obs['dpt_pseudotime'], df=ns_df)
            if len(null_obs) > 1:
                null_obs = null_obs[null_obs != 'dpt_pseudotime']
                Nnullobs = len(null_obs)
                null_obs = adata.obs[null_obs].to_numpy()
                null_obs = null_obs.astype(float)
                null_obs = null_obs.reshape(Ncells, len(Nnullobs))
                null_obs = np.concatenate((null_obs, pt_design), axis=1)
            else:
                null_obs = pt_design
        else:
            null_obs = adata.obs[null_obs].to_numpy()
    null_obs = null_obs.astype(float)
    
    pval_all = np.array([])
    for col_idx in range(Ngenes):
        pval = likelihood_ratio_test(X_alt=alt_obs, y=yall[:,col_idx], X_null=null_obs)
        pval_all = np.append(pval_all, pval)
    pval_all[np.isnan(pval_all)] = 1
    _, padj, _, _ = multi.multipletests(pval_all, method='fdr_bh')
    results = pd.DataFrame(data=np.concatenate((pval_all.reshape(Ngenes,1), padj.reshape(Ngenes,1)), axis=1), index=GN, columns=['pval', 'padj'])
    results = results.iloc[np.argsort(results['padj']),:]
    return results

def smoothing_fun(x, w=0.1, sigma=20):
    '''
    Apply Guassian window for smoothing
    
    Input
    -----
    x: an array of values
    w: percentage of length of x as the number of points in the output window, M=w*len(x), default: 0.1
    sigma: the standard deviation, default: 20
    
    Returns
    -----
    Smoothed values
    
    '''
        
    window = gaussian(np.floor(w*len(x)), sigma)
    smoothed = np.convolve(x, window / window.sum(), mode='same')
    return smoothed

def plot_genes_along_pt(adata, genes, pt_obs='dpt_pseudotime', figsize=(6,4), smooth=True, save=None, **kwargs):
    '''
    General plotting function for genes along pseudotime
    
    Input
    -----
    adata: adata obj
    genes: genes to plot
    pt_obs: PT .obs key
    figsize: size of figure, default: (6,4)
    smooth: if smooth using gaussian window, default: True
    save: name for saving the plot, default: None
    **kwargs: other parameters in the smoothing function
    
    Returns
    -----
    A general plot with (smoothed) gene expression pattern along PT
    
    '''
        
    pt = adata.obs[pt_obs].values
    pt_idx = np.argsort(pt)
    pt = pt[pt_idx]
    difG = np.setdiff1d(genes, adata.raw.var_names)
    if len(difG) != 0:
        raise ValueError(f'Genes: {difG} do not exist in adata.raw.')
    
    plt.figure(figsize=figsize)
    for g in genes: 
        gExp = adata.raw[pt_idx, g].X.flatten()
        if smooth:
            gExp = smoothing_fun(gExp, **kwargs)
        gExp = (gExp-np.min(gExp))/(np.max(gExp)-np.min(gExp))
        plt.plot(pt, gExp, label=g)
    plt.grid(False)
    plt.xlabel('diffusion pseudotime')
    plt.yticks([0,1], ('min', 'max'))
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    
    if save is not None:
        plt.savefig(save)
        
######## BELOW is for projection #########################
## quick neighbours detection for dpt PT calculation ##
from scipy.sparse import coo_matrix
from umap.umap_ import fuzzy_simplicial_set
from sklearn.metrics import pairwise_distances
def _get_sparse_matrix_from_indices_distances_umap(knn_indices, knn_dists, n_obs, n_neighbors):
    '''
    This code is copied from scanpy!!!
    '''
    
    rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    vals = np.zeros((n_obs * n_neighbors), dtype=np.float64)

    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] == -1:
                continue  # We didn't get the full knn for i
            if knn_indices[i, j] == i:
                val = 0.0
            else:
                val = knn_dists[i, j]

            rows[i * n_neighbors + j] = i
            cols[i * n_neighbors + j] = knn_indices[i, j]
            vals[i * n_neighbors + j] = val

    result = coo_matrix((vals, (rows, cols)),
                                      shape=(n_obs, n_obs))
    result.eliminate_zeros()
    return result.tocsr()

def compute_connectivities_umap(
    knn_indices, knn_dists,
    n_obs, n_neighbors, set_op_mix_ratio=1.0,
    local_connectivity=1.0,
):

    """\
    This code is copied from scanpy!!!
    
    This is from umap.fuzzy_simplicial_set [McInnes18]_.
    Given a set of data X, a neighborhood size, and a measure of distance
    compute the fuzzy simplicial set (here represented as a fuzzy graph in
    the form of a sparse matrix) associated to the data. This is done by
    locally approximating geodesic distance at each point, creating a fuzzy
    simplicial set for each such point, and then combining all the local
    fuzzy simplicial sets into a global one via a fuzzy union.
    """

    X = coo_matrix(([], ([], [])), shape=(n_obs, 1))
    connectivities = fuzzy_simplicial_set(
        X,
        n_neighbors,
        None,
        None,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
        set_op_mix_ratio=set_op_mix_ratio,
        local_connectivity=local_connectivity,
    )

    if isinstance(connectivities, tuple):
        # In umap-learn 0.4, this returns (result, sigmas, rhos)
        connectivities = connectivities[0]

    distances = _get_sparse_matrix_from_indices_distances_umap(
        knn_indices, knn_dists, n_obs, n_neighbors
    )

    return distances, connectivities.tocsr()


def quick_neighbors(comb, metric='euclidean', n_neighbors = 10, random_state = 0):
    '''
    A quick neighbours calculation between reference data and the new data.
    If reference data/new data has batch effect, it needs to be firstly
    corrected using fastMNN, which does correction based on the PCA space.
    The corrected PCA will be used for neighbours calculation.
    
    Assumptions:
    1) Cells from ref data will only be connected to each other
    2) Cells from new data will only be connected to Cells from ref data
    
    This function generates the same output as sc.pp.neighbors function in scanpy.
    
    Input
    -----
    comb: combined adata obj. done by adata_ref.concatenate(adata_new), adata_ref must always be at the front
    metric: method for distance calculation, for details, please see scanpy sc.pp.neighbors function, default: 'euclidean'
    n_neighbors: number of neighbours to consider, default: 10
    random_state: random seed, default: 0
    
    Returns
    -----
    comb.uns['neighbors'] same as outputs from sc.pp.neighbors
    
    '''
    
    if 'X_pca' not in comb.obsm_keys():
        raise ValueError('PCA needs to be calculated and combined first')
    X = comb.obsm['X_pca']
    D = pairwise_distances(X, metric=metric)
    nCell_total = X.shape[0]
    nCell_ref = sum(comb.obs['batch'] == '0')
    nCell_new = sum(comb.obs['batch'] == '1')
    
    sample_range = np.arange(nCell_total)[:, None]
    knn_indices = np.argpartition(D[:,0:nCell_ref], n_neighbors-1, axis=1)[:, :n_neighbors]
    knn_indices = knn_indices[sample_range, np.argsort(D[sample_range, knn_indices])]
    knn_dists = D[sample_range, knn_indices]
    distances, connectivities=compute_connectivities_umap(knn_indices, knn_dists, nCell_total, n_neighbors)

    comb.uns['neighbors'] = {}
    comb.uns['neighbors']['params'] = {'n_neighbors': n_neighbors, 'method': 'umap'}
    comb.uns['neighbors']['connectivities'] = connectivities
    comb.uns['neighbors']['distances'] = distances
    
## quick calculation of umap and umap projection ######
import umap
def quick_umap(adata_ref, n_neighbors=10, min_dist: float = 0.5, spread: float = 1.0,
              n_components: int = 2, alpha: float = 1.0, a = None, b=None,
               negative_sample_rate: int = 5, init_coords = 'spectral', 
               random_state = 0, **kwargs
              ):
    
    '''
    Calculate a UMAP for the reference data using the umap python package.
    
    Not using sc.tl.umap is due to umap model needs to be saved first as output.
    
    Most of the default parameters are from scanpy sc.tl.umap function.
    
    Input
    -----
    adata_ref: reference anndata object, needs to have 'X_pca' in .obsm and 'neighbors' in .uns
    n_neighbors: Number of neighbors to be considered, default: 10
    min_dist: The effective minimum distance between embedded points, default: 0.5
    spread: The effective scale of embedded points. In combination with `min_dist`
        this determines how clustered/clumped the embedded points are. Default: 1.0
    n_components: The number of dimensions of the embedding, default: 2 
    alpha: The initial learning rate for the embedding optimization, default: 1.0
    a: More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`. Default: None
    b: More specific parameters controlling the embedding. If `None` these
        values are set automatically as determined by `min_dist` and
        `spread`. Default: None
    negative_sample_rate: The number of negative edge/1-simplex samples to use per positive
        edge/1-simplex sample in optimizing the low dimensional embedding. Default: 5
    int_coords: How to initialize the low dimensional embedding. Called `init` in the
        original UMAP. Default: 'spectral'
    random_stat: random seed, default: 0
    **kwargs: other parameters taken in the umap.umap_.UMAP function
    
    Returns
    -----
    A umap model
    The calculate umap coordicates for ref data is saved in data_ref.obsm['X_umap']
    
    '''
    if 'X_pca' not in adata_ref.obsm_keys():
        raise ValueError('Need to calculate PCA first')
    
    if 'neighbors' not in adata_ref.uns_keys():
        raise ValueError('Need to calcualte neighbors first')
    
    if a is None or b is None:
        a, b = umap.umap_.find_ab_params(spread, min_dist)
    neighbors = adata_ref.uns['neighbors']
    neigh_params = neighbors['params']
    umap_ref = umap.umap_.UMAP(
        n_neighbors = n_neighbors,
        n_components = n_components,
        learning_rate = alpha,
        a = a,
        b = b,
        negative_sample_rate = negative_sample_rate,
        init = init_coords,
        random_state = random_state,
        output_metric = neigh_params.get('metric', 'euclidean'),
        output_metric_kwds = neigh_params.get('metric_kwds', {}),
        **kwargs
    )
    X = adata_ref.obsm['X_pca']
    X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
    X_umap_fit = umap_ref.fit(X_contiguous)
    X_umap = X_umap_fit.embedding_
    adata_ref.obsm['X_umap'] = X_umap
    return umap_ref

def quick_umap_proj(adata_new, umap_ref):
    '''
    Project new data onto the existing reference data umap based on pca space
    
    Input
    -----
    adata_new: new anndata object, needs to have 'X_pca' in .obsm
    umap_ref: umap model as the output of the 'quick_umap' function
    
    Returns
    -----
    Projected umap coordicates of the new data
    
    '''
    
    if 'X_pca' not in adata_new.obsm_keys():
        raise ValueError('Need to calculate PCA first')
    
    X1 = adata_new.obsm['X_pca']
    X1_contiguous = np.ascontiguousarray(X1, dtype=np.float32)
    X1_umap = umap_ref.transform(X1_contiguous)
    return X1_umap

################### BELOW is for plotting 3d ###################################
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.1, 0.5, 0.5),
                 (0.2, 0.0, 0.0),
                 (0.4, 0.2, 0.2),
                 (0.6, 0.0, 0.0),
                 (0.8, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),
        'green':((0.0, 0.0, 0.0),
                 (0.1, 0.0, 0.0),
                 (0.2, 0.0, 0.0),
                 (0.4, 1.0, 1.0),
                 (0.6, 1.0, 1.0),
                 (0.8, 1.0, 1.0),
                 (1.0, 0.0, 0.0)),
        'blue': ((0.0, 0.0, 0.0),
                 (0.1, 0.5, 0.5),
                 (0.2, 1.0, 1.0),
                 (0.4, 1.0, 1.0),
                 (0.6, 0.0, 0.0),
                 (0.8, 0.0, 0.0),
                 (1.0, 0.0, 0.0))}

def plot_3d(adata_ref, obs_key, adata_new=None, obsm_key='X_diffmap', ncols=4,figsize=(6,6), 
            alpha=0.5, azim=250,elev=30, markersize=1,components=[1,2,3], cmap=None, save=None):
    '''
    This function is temporarily used for 3d plot as scanpy projection='3d' function does not
    work properly with the latest version due to matplotlib problem.
    
    Input
    -----
    adata_ref: reference anndata object
    obs_key: the obs key to plot, if adata_new is None, then it is the obs key from adata_ref, otherwise, it is from adata_new
    adata_new: new anndata object, if this is not None, then adata_ref will be plotted as background
    obsm_key: obsm layout, default: X_diffmap
    ncols: number of columns for each row, default: 4
    figsize: figure size, default: (6,6)
    alpha: transparency, default: 0.5
    azim: rotation parameter, default: 250
    elev: rotation parameter, default: 30
    markersize: point size, default: 1
    components: components to plot, default: 1,2,3
    cmap: color map, if None, using default rainbow color
    save: figure name for saving, default: None
    
    Returns
    -----
    3d scatter plot showing the layouts
    
    '''
    nkey = len(obs_key)
    if nkey <=4:
        ncols = nkey
    nrows = int(np.ceil((nkey)/ncols))
    fig = plt.figure(figsize=figsize)
    for nk in range(nkey):
        k = obs_key[nk]
        ax = fig.add_subplot(nrows, ncols, nk+1, projection='3d')
        ax.view_init(azim=azim, elev=elev)
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        
        if adata_new is None:
            adata = adata_ref.copy()
        else:
            adata = adata_new.copy()
            dm_ref = adata_ref.obsm[obsm_key]
            ax.plot(dm_ref[:,components[0]],dm_ref[:,components[1]],dm_ref[:,components[2]], '.', markersize=markersize, c='#d3d3d3', label = 'Ref_data', alpha =alpha)

        dm_new = adata.obsm[obsm_key]
        if k in adata.obs_keys():
            obs_term = np.array(adata.obs[k].values)
        elif k in adata.raw.var_names:
            obs_term = adata.raw[:, k].X.flatten()
        if obs_term.dtype == 'float' or obs_term.dtype == int or obs_term.dtype == 'float32':
            if cmap is None:
                cmap = LinearSegmentedColormap('my_colormap',cdict,256) 
            #print(np.array(obs_term.values))
            conti_fig = ax.scatter(dm_new[:,components[0]],dm_new[:,components[1]],dm_new[:,components[2]], '.', s=markersize, c=obs_term, cmap=cmap, alpha =alpha)
            fig.colorbar(conti_fig, shrink=0.5)
            ax.set_title(k)
        
        else:
            obs_term = adata.obs[k].astype('category')
            cats = obs_term.cat.categories
            if k+'_colors' in adata.uns_keys():
                color_pal = adata.uns[k+'_colors']
            else:
                color_pal = sc.pl.palettes.default_20[0:len(cats)]
            for i in range(len(cats)):
                #print(cats[i])
                idx = obs_term==cats[i]
                ax.plot(dm_new[idx,components[0]],dm_new[idx,components[1]],dm_new[idx,components[2]], '.', markersize=markersize, c=color_pal[i], label = cats[i], alpha =alpha)
            ax.set_title(k)
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
    plt.tight_layout()
    if save is not None:
        plt.savefig(save)
        

###################### BELOW is for pathway analysis ############################################
def pathway_score_cal(adata, DBcont, minGS=5, maxGS=500):
    '''
    This function is to calculate geometric mean for each terms in the database for each cell.
    
    Input
    -----
    adata: the anndata object
    DBcont: the gmt read in file, can be either downloaded from GSEA website or user defined
    minGS: min number of genes in the gene set, default: 5
    maxGS: max number of genes in the gene set, default: 500
    
    Returns
    -----
    Exp pd DataFrame for database, rows are cells and columns are DB terms.
    
    '''
    expArray = []
    pnames = np.array([])
    TotalGenes = [x.upper() for x in adata.raw.var_names]
    for l in range(len(DBcont)):
        lcont = DBcont[l].split('\t')
        pathway = lcont[0]
        DBGenes = [x.upper() for x in lcont[2:]]
        DBGenes = np.in1d(TotalGenes, DBGenes)
        if ((np.sum(DBGenes) < minGS) or (np.sum(DBGenes) > maxGS)): 
            continue
        else:
            Exp = np.mean(adata.raw[:, DBGenes].X, axis=1)
            expArray.append(Exp)
            pnames = np.append(pnames, pathway)
    expArray=pd.DataFrame(np.vstack(expArray).T)
    expArray.index = adata.obs_names
    expArray.columns = pnames
    return expArray

import scipy.stats as stats
import statsmodels.stats.multitest as multi
def pathway_analysis(GL, TotalGenes, DBcont, minGS=5, maxGS=500):
    '''
    This function is to calculate the enriched database terms for a given gene set using hypergeometric test.
    
    Input
    -----
    GL: a given gene set
    TotalGenes: total number of genes, can be found from len(adata.raw.var_names)
    DBcont: the gmt read in file, can be either downloaded from GSEA website or user defined
    minGS: min number of genes in the gene set, default: 5
    maxGS: max number of genes in the gene set, default: 500
    
    Returns
    -----
    pd DataFrame for each term of the databse, including:
    Pathway: name of the pathway
    k: number of overlapped genes between the gene list and the pathway
    M: number of total genes
    n: numbber of genes in the pathway
    N: number of genes in the gene list
    pval: hypergeometric p value
    OLGenes: overlapped genes
    padj: adjusted p value by BH method
    
    '''
    TotalGenes = [x.upper() for x in TotalGenes]
    M = len(TotalGenes)
    GL = [x.upper() for x in GL]
    
    saveDict = {}
    saveDict['Pathway'] = []
    saveDict['k'] = []
    saveDict['M'] = []
    saveDict['n'] = []
    saveDict['N'] = []
    saveDict['pval'] = []
    saveDict['OLGenes'] = []
    for l in range(len(DBcont)):
        lcont = DBcont[l].split('\t')
        pathway = lcont[0]
        DBGenes = [x.upper() for x in lcont[2:]]
        DBGenes = np.intersect1d(DBGenes, TotalGenes)
        n = len(DBGenes)
        N = len(GL)
        OLGenes = np.intersect1d(DBGenes, GL)
        k = len(OLGenes)
        if ((n >=minGS) or (n <= maxGS)):
            #print('k:'+str(k-1)+' M:'+str(M)+' n:'+str(n)+' N:'+str(N))
            pval = stats.hypergeom.sf(k-1, M, n, N)
            saveDict['Pathway'].append(pathway)
            saveDict['k'].append(k)
            saveDict['M'].append(M)
            saveDict['n'].append(n)
            saveDict['N'].append(N)
            saveDict['pval'].append(pval)
            saveDict['OLGenes'].append(';'.join(OLGenes))
    pA = np.array(saveDict['pval'])
    pA[np.isnan(pA)] = 1
    _, padj, _, _ = multi.multipletests(pA, method='fdr_bh')
    saveDict['padj'] = padj
    saveDF = pd.DataFrame.from_dict(saveDict)
    saveDF = saveDF.sort_values(by='padj')
    return saveDF

