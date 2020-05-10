#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
import anndata
from scipy.stats import chi2
from scipy.sparse import issparse
import statsmodels.stats.multitest as multi
import statsmodels.api as sm

try:
    from scanpy import logging as logg
except ImportError:
    pass


def read_in_files(Indir, ftable_loc, method = 'FeatureCount'):
    '''
    Read in STAR aligned files. 
    
    Input
    -----
    Indir: FeatureCount output folder path or HTSeqcount output file path
    ftable_loc: Gene feature table path, ftable must have two columns with ['Gene Name', 'Gene Type']
    method: Counting method, can be either 'FeatureCount' or 'HTSeqcount', default: FeatureCount
    
    Returns
    -----
    An anndata object
    
    '''
    
    if method == 'FeatureCount':
    # first tidy up X
        X = pd.read_csv(Indir+'/fcounts/fcounts.txt', delimiter='\t', index_col=0, skiprows=1)
        GT = X.iloc[:-92,0:5].copy()
        X = X.iloc[:,5:].transpose().copy()
        X['obs_names'] = [re.search('(SLX-\d+\.\w\d+_\w\d+)', x).group(0) for x in X.index]
        X = X.groupby('obs_names').sum()
        ERCC = X.iloc[:,-92:].copy()
        X = X.iloc[:,:-92].copy()

        # Then work on the QC
        QC = pd.read_csv(Indir+'/fcounts/fcounts.txt.summary', delimiter='\t', index_col=0).transpose()
        QC['obs_names'] = [re.search('(SLX-\d+\.\w\d+_\w\d+)', x).group(0) for x in QC.index]
        QC = QC.groupby('obs_names').sum()
        QC = QC.iloc[:,(np.sum(QC, axis=0)!=0).values].copy()
        QC.columns = ['QC_'+x for x in QC.columns]
        QC = QC[[x for x in QC.columns if 'Unassigned' in x]]
    elif method == 'HTSeqcount':
        X = pd.read_csv(Indir, delimiter='\t', index_col=0).transpose()
        ERCC = X.loc[:,['ERCC' in x for x in X.columns]].copy()
        QC = X.loc[:,['__' in x for x in X.columns]].copy()
        QC.columns = [x.replace('__','') for x in QC.columns]
        QC.columns = ['QC_'+x for x in QC.columns]
    else:
        raise ValueError('method can only be either FeatureCount or HTSeqcount.')
    
    # Combine feature table
    ftable = pd.read_csv(ftable_loc, delimiter='\t', index_col=0, header=None)
    ftable.columns = ['Gene Name', 'Gene Type']
    if method == 'FeatureCount':
        FeatureTable = pd.concat([ftable,GT], axis=1)
    else:
        FeatureTable = ftable
    
    # New write in anndata frame
    adata = anndata.AnnData(X=X, var=FeatureTable, obs=QC)
    adata.obsm['ERCC'] = ERCC
    adata.var['Ensembl_ID'] = adata.var_names
    adata.var_names = adata.var['Gene Name']
    adata.var_names_make_unique()
    return adata

def smartseq_qc(adata, cutoff=[np.log10(2*(10**5)), 0, 0.2, 4000, 0.2, 0.2, float('nan'), float('nan'), float('nan'), float('nan')],
            MTpattern = 'mt-', ncols=4, figsize=(10,7), s=10, title=None, save=None):
    '''
    Do a bglab equivalent quality control.
    
    Input
    -----
    adata: anndata object
    cutoff: QC cutoffs, 6 in total, ordered by ['nMapped (log10)', 'nNuclear (log10)', 'fGenes:nTotal', 'nHCGenes', 'mito:nGenes', 'nERCC:nMapped']
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
    
    mito_genes = [name for name in adata.var_names if name.startswith(MTpattern)]
    mitoCNT = np.sum(adata[:,mito_genes].X, axis=1).copy()
    nuclearCNT = np.sum(adata[:,~np.in1d(adata.var_names,mito_genes)].X, axis=1).copy()
    erccCNT = np.sum(adata.obsm['ERCC'], axis=1).values
    qcNames = [x for x in adata.obs_keys() if 'QC' in x]
    qcCNT = np.sum(adata.obs[qcNames], axis=1).values
    nTotal = mitoCNT + nuclearCNT + erccCNT + qcCNT
    nMapped = mitoCNT + nuclearCNT + erccCNT
    nGenes = mitoCNT + nuclearCNT
    nHCGenes = np.sum(adata[:,~np.in1d(adata.var_names,mito_genes)].X.T*(10**6)/nuclearCNT > 10, axis=0)
    QCdata = {}
    QCdata['nMapped (log10)'] = np.log10(nMapped)
    QCdata['nNuclear (log10)'] = np.log10(nuclearCNT)
    QCdata['fGenes:nTotal'] = nGenes/nTotal
    QCdata['nHCGenes'] = nHCGenes
    QCdata['mito:nGenes'] = mitoCNT/nGenes
    QCdata['nERCC:nMapped'] = erccCNT/nMapped
    QCdata['nNuclear:nMapped'] = nuclearCNT/nMapped
    for qcIndex in qcNames:
        QCdata[qcIndex.replace('_Unassigned', '')+':nTotal'] = adata.obs[qcIndex]/nTotal
    
    # cells failed QC
    compara = ['<','<','<','<','>','>','<','>','>','>']
    failed = []
    for i in range(len(QCdata.keys())):
        if not np.isnan(cutoff[i]):
            if compara[i] == '<':
                failed.append(QCdata[list(QCdata.keys())[i]] < cutoff[i])
            else:
                failed.append(QCdata[list(QCdata.keys())[i]] > cutoff[i])
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
        if not np.isnan(cutoff[i]):
            ax[rowidx, colidx].axhline(y=cutoff[i], color='orange', linestyle='dashed')
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
    adata.obs['percent_mito'] = mitoCNT/nGenes
    adata.obs['n_genes'] = np.sum(adata.X > 0, axis=1)
    
    return adata[~failed_idx,:].copy()

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
    #print(meta.columns)
    meta.columns = ['Gottgens_ID', 'Sequencing_identifier', 'Plate_number',
       'Position_in_96_well_plate_sorted', 'Position_in_96_well_plate_RNAseq',
       'FACs', 'Unique_sample_descriptor', 'Details', 'Cell_type_general',
       'Cell_type_subtype', 'Owner', 'Species', 'Sequencing_index',
       'Sequencing_facility_index', 'Average_pooled_library_length_bp', 'Pool_size',
       'Number_of_lanes']
    #print(meta.columns)
    return(meta)

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

def normalise_data(adata, reCalSF=True, copy=False):
    '''
    Normalisation using size factors estimated by DESeq2 method
    
    Input
    -----
    adata: an anndata object
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
        sf_genes = est_size_factor(adata.X)
        adata.obs['sf_gene'] = sf_genes
        if 'ERCC' in adata.obsm_keys():
            print('Calculate SF for erccs:')
            sf_ercc = est_size_factor(adata.obsm['ERCC'])
            adata.obs['sf_ercc'] = sf_ercc
    else:
        if 'sf_gene' not in adata.obs_keys():
            raise ValueError('sf_gene is not found in .obs, please set reCalSF=True.')
    
    adata.X = np.log1p(adata.X/adata.obs['sf_gene'][:,None])
    if 'ERCC' in adata.obsm_keys():
        adata.obsm['ERCC_norm'] = np.log1p(adata.obsm['ERCC']/adata.obs['sf_ercc'][:,None])
    if copy:
        return adata.copy()
    

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
    link_func = sm.genmod.families.links.identity
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
    ax.set_yscale('log',basey=10)
    ax.set_xscale('log',basex=10)
    ax.set_ylabel(r'$\sigma^{2}/\mu^{2}$')
    ax.set_xlabel(r'$\mu$')
    ax.grid(False)

    plt.scatter(g_df.loc[g_df['highVar'],:]['mean'], g_df.loc[g_df['highVar'],:]['cv2'], color='red', s=s)
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    xpred = np.linspace(start=ax.get_xlim()[0], stop=ax.get_xlim()[1], num=1000)
    if useERCC:
        plt.scatter(e_df['mean'], e_df['cv2'], s=s, color='blue')
        y1 = a1t/xpred + a0
        y2 = psi/xpred + a0 + minBiolDisp
        plt.plot(xpred[y1<=ylim[1]], y1[y1<=ylim[1]], color='red')
        plt.plot(xpred[y1<=ylim[1]], y2[y1<=ylim[1]], color='red', linestyle='dashed')
    else:
        y1 = psi/xpred + a0
        plt.plot(xpred[y1<=ylim[1]], y1[y1<=ylim[1]], color='red', linestyle='dashed')
    
    if save is not None:
        plt.savefig(save)
        
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
        adata_sub = adata_sub.todense()
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

