#!/home/couturier/.local/share/canopy/edm/envs/User/bin/python
Help=\
"""
process_gbm2.py
Author: Charles Couturier
Created: 2017/3/27
Last update: 2019/1/29

This program is made to extract appropriate signal from gene-barcode matrices in
10X format. Running from command line, it creates a
.mat file and a .h5 file of the resulting gbm as well as analysis plots. 
Standard output can be used to create a report of the processes. 

Steps are as follows:
-Loading gbm
-Removing bad cells (cells with high MT genes proportions or low gene count)
-Calculating and adjusting for size factor
-Extracting signal based on Fano factors
-Normalizing (log and z-score)
-Running pca for visualization

Parameters to specify:

'-f': input file
'-d': input directory
'-o': output directory
'-s': size factor method ( 'median' (DESeq) OR 'total' (TPM) )
'-t': thresholds (see below)
'-g': genome used for mapping by cellranger (eg. GRCh38)
'-p': plot (True or False - case insensitive)
'-m': comma seperated list of markers (symbols) to plot in sample PCA plots
'-z': file containing barcodes for outliers (separate barcodes with newline)
'-i': comma separated list of genes for which to force inclusion (if >0 count)
'-h': print this help text

Only input file OR directory need be specified. If both are specified,
input file is used. If not directory is used (10X outs directory); looking for
h5 file first; if not found, looking for mtx; if not found, looking for csv
last.  

Thresholds, as specified by '-t', are given by a comma seperated list of values.
Not all have to be specified. Add commas for all values to skip (eg. for second
threshold: ,3000; for first and last: 0.12,,,,1.5)  
In order, they are
1) mitochondrial proportion cutoff for acceptable cell (typically 0.12 to 0.3)
2) genes per cell cutoff (typically 2000-3000)
3) UMI per cell cutoff (10X sets it at 2000)
4) cutoff for expression ratio of highly expressed genes to exclude from 
size factor determination
5) minimum number of cells which must express gene in order to be considered 
by Fano plot
6) Fano factor threshold 
Default: 0.12,2000,0,0.01,10,1.5 

2019/1/29: update to read latest cellranger output

Example of syntax:
process_gbm2.py -d outs -o test -m SOX2,GFAP,OLIG2 > test.log   

"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import scipy as sp
import scipy.io as sio
import math
import tables
import csv
import copy
import warnings

#--------------------------------------------------------------------------
# Miscellaneous functions and class

class GeneBCMatrix:
    
    def __init__(self,ensembl, ensembl2symb, symb2ensembl, barcodes, matrix):
        self.ensembl = ensembl                      # list of all genes by ensembl
        self.e2s = ensembl2symb                     # dict ensembl to symb
        self.s2e = symb2ensembl                     # dict symb to ensembl
        self.barcodes = barcodes                    # cell barcodes
        self.matrix = matrix.astype(float)          # matrix with data
        self.raw_counts = None                      # raw data, none at first
        self.log = None                             # log data, none at first
        self.bd = {b:i for i,b in enumerate(self.barcodes)}# barcode dictionary
        if len(self.bd)==len(self.barcodes):
            print 'All barcodes are unique'
        else:
            print 'Duplicate barcodes exist'
        
    def __getitem__(self,slices):
        if type(slices)==tuple:
            if len(slices)>2:
                sys.stderr.write\
                    ('gbm has two dimensions - more than two given as arg\n')
                pass
            elif len(slices)==2:
                self.matrix = self.matrix[slices[0],slices[1]]
                self.barcodes = self.barcodes[slices[0]]
                self.ensembl = self.ensembl[slices[1]]
        else: 
            self.matrix = self.matrix[slices]
            self.barcodes = self.barcodes[slices]
        self.bd = {b:i for i,b in enumerate(self.barcodes)} # update bd
            

# inverse dictionary from ensembl --> symbol, to symbol --> ensembl 
def inv_dict(my_map):
    inv_map = {}
    not_unique = []
    for k, v in my_map.iteritems():
        inv_map[v] = inv_map.get(v, [])
        if inv_map[v] != []:
            not_unique.append(v)
        inv_map[v].append(k)
    if not_unique != []:
        print 'No unique mapping for following symbols:'
        print not_unique
    #else:
     #   del inv_map
      #  inv_map = {v: k for k, v in my_map.iteritems()}
    return inv_map

# get name of gbm file based on availability
def get_gbm_file(dir_, input_file, genome):
    gbm_type = 'filtered_feature_bc_matrix'
    h5_file = os.path.join(dir_, gbm_type+'.h5')
    mtx_file = os.path.join(dir_,gbm_type, genome,'matrix.mtx')
    
    if os.path.isfile(input_file) == False:
        if os.path.isfile(h5_file):
            input_file = h5_file
        elif os.path.isfile(mtx_file):
            input_file = mtx_file
        else:
            input_file = ''
    
    return input_file

# get gbm from h5 file - adapted from script provided by 10X
def get_matrix_from_h5(filename, genome):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')

	feature_group = f.get_node(mat_group, 'features')
        ensembl = getattr(feature_group, 'id').read()
        symb = getattr(feature_group, 'name').read()

	barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()

        matrix = sp.sparse.csc_matrix((data, indices, indptr), shape=shape)
        matrix = matrix.T
        return ensembl, symb, barcodes, matrix


# ------------------------------------------------------------------------------
# load gbm csv file. 'gbm_file' is a file directory    
def load_gbm(gbm_file, genome):
    sys.stderr.write('Loading gbm... ')
    
    if gbm_file == '':
        raise TypeError('gbm_file not found')
        
    if gbm_file.endswith('.h5'):
        ensembl, symb, barcodes, matrix = get_matrix_from_h5(gbm_file, genome)
        
    if gbm_file.endswith('.mtx'):
        human_matrix_dir = os.path.dirname(os.path.realpath(gbm_file))
        matrix = sp.io.mmread(os.path.join(human_matrix_dir, "matrix.mtx"))
        matrix = sp.sparse.csr_matrix(matrix.T)
        
        genes_path = os.path.join(human_matrix_dir, "genes.tsv")
        ensembl = np.array(\
            [row[0] for row in csv.reader(open(genes_path), delimiter="\t")] )
        symb = np.array(\
            [row[1] for row in csv.reader(open(genes_path), delimiter="\t")] )
 
        barcodes_path = os.path.join(human_matrix_dir, "barcodes.tsv")
        barcodes = np.array(\
            [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")])
        
    if gbm_file.endswith('.csv'):
        gbm = pd.read_csv(gbm_file,header=None,sep=',',index_col=0).transpose()
        symb = gbm.iloc[0]
        ensembl = gbm.columns
        gbm = gbm.drop(gbm.index[0])
        matrix = gbm.as_matrix().astype(float)
        matrix = sp.sparse.csc_matrix(matrix)
       
    ensembl2symb = {ensembl[i]:symb[i] \
        for i in range(len(ensembl))}
    symb2ensembl = inv_dict(ensembl2symb) # outputs a dict of lists   
    sys.stderr.write('done\n')
    return GeneBCMatrix(ensembl, ensembl2symb, symb2ensembl, barcodes, matrix)

    
#------------------------------------------------------------------------------- 
# removing artifactual cells, based on MT gene expression and total genes

def find_badcells(gbm,MT_cutoff=0.12, gpc_cutoff=2000, UMIpc_cutoff=0,\
    do_plot='True'):
    
    sys.stderr.write('Removing bad cells... ')
    
    # remove cells with high MT genes transcription 
    MT_symb = ['MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6',\
        'MT-CYB','MT-CO1','MT-CO2','MT-CO3','MT-ATP6','MT-ATP8','MT-RNR2']
    MT = []
    for mt in MT_symb:
        try:
            MT.extend(gbm.s2e[mt])
        except KeyError:
            pass
            
    find = lambda searchList, elem: [np.where(searchList==e)[0][0] for e in elem]
    pMT = gbm.matrix[:,find(gbm.ensembl,MT)].sum(axis=1)/gbm.matrix.sum(axis=1)
    pMT = np.squeeze(np.asarray(pMT.T))
    plotpMT = np.sort(pMT)[::-1]
    
    # save plot for MT distribution
    if do_plot.lower()=='True'.lower():
        plt.figure()
        plt.plot(plotpMT)
        try:
            plt.axvline(x=np.where(plotpMT>MT_cutoff)[0][-1],color='k')
        except IndexError:
            pass
        plt.ylim(0,1)
        plt.xlabel('cell index')
        plt.ylabel('MT proportion')
        plt.savefig('MT_distribution.eps', format='eps')
    
    # genes_per_cell
    genes_per_cells = np.squeeze(np.asarray(\
        gbm.matrix.astype(bool).sum(axis=1)))
    pltgpc = np.sort(genes_per_cells)[::-1]
    
    plt.figure()
    plt.plot(pltgpc,c='blue')
    try:
        plt.axvline(x=np.where(pltgpc>gpc_cutoff)[0][-1],color='k')
    except IndexError:
        pass
    plt.title('genes per cell')
    plt.xlabel('cell')
    plt.ylabel('number of different genes')
    plt.savefig('genes_per_cell.eps',bbox_inches='tight',\
        format='eps',dpi=100)  
    
    # UMI per cell  
    UMI_per_cell = np.squeeze(np.asarray(gbm.matrix.sum(axis=1)))
    pltUpc = np.sort(UMI_per_cell)[::-1]
    plt.figure()
    plt.subplot(2,1,1)
    plt.hist(UMI_per_cell, bins=100)
    plt.ylabel('number of cells')
    plt.xlabel('UMI count')
    plt.subplot(2,1,2)
    plt.plot(pltUpc)
    try:
        plt.axvline(x=np.where(pltUpc>UMIpc_cutoff)[0][-1],color='k')
    except IndexError:
        pass
    plt.ylabel('UMI count')
    plt.xlabel('cells')
    plt.savefig('UMI_per_cell.eps',bbox_inches='tight',\
        format='eps',dpi=100)
            
    sys.stderr.write('done\n')
    return (pMT<MT_cutoff) & (genes_per_cells>gpc_cutoff) & \
        (UMI_per_cell>UMIpc_cutoff)

# ------------------------------------------------------------------------------
# find outliers, return cells to keep
def remove_outliers(gbm,outliers_file):
    sys.stderr.write('Removing outliers...')
    keep_cells = np.ones((gbm.matrix.shape[0]),dtype=bool)
    #if outliers_file=='':
    if os.path.isfile(outliers_file)==False:
        sys.stderr.write('none given\n')
    else:
        with open(outliers_file) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=' ')
            outliers = [row[0] for row in spamreader if row]
        print 'Removed cells listed in %s'%outliers_file
        print '%i outliers to remove'%len(outliers)
        olist=[]
        for o in outliers:
            try:
                olist.append(gbm.bd[o])
            except KeyError:
                print '%s not found'%o
        index = np.array(olist)
        print 'Removed %i outliers'%len(olist)
        #index = np.array([gbm.bd[o] for o in outliers])
        if olist!=[]: keep_cells[index] = False    
        sys.stderr.write('done\n')
    
    return keep_cells
            
# ------------------------------------------------------------------------------
# normalize each cell by size. Methods: median with size factor (DESeq), 
#   total (TPM)
def size_factor(gbm,method='total',thresh=0.01):
    sys.stderr.write('Calculating size factor...')
    data = gbm.matrix
    if method=='median':
        print ' (median UMI)... ',
        s = ( data+1/(sp.stats.mstats.gmean(data+1,axis=0)) ).median(axis=1)
        
    elif method=='total':
        print 'Using total UMI for size factor'
        data = data.astype(float)
        ER = (data.T/data.sum(axis=1).T).T
        tdata=data
        tdata[ER>thresh]=0
        #mean_tUMI = tdata.sum(axis=1).mean()
        mean_tUMI = 1e5
        s = tdata.sum(axis=1)/mean_tUMI
        
    else:
        raise TypeError('unknown method')
    
    sys.stderr.write('done\n')
    return s
    
#-------------------------------------------------------------------------------
# Find highly variable genes. Uses Fano factor.

def get_signal(gbm, cthresh=10, fthresh=1.5, do_plot='True'):
    sys.stderr.write('Extracting signal... ')
    data = gbm.matrix.todense()
    
    # remove genes with less than 'cthresh' cells with at least one count
    cell_count = (data>=1).sum(axis=0)
    cell_count = np.squeeze(np.asarray(cell_count))
    data = data[:,cell_count>=cthresh]
    
    # calculate Fano distribution
    gene_var = np.squeeze(np.asarray(data.var(axis=0)))
    gene_mean = np.squeeze(np.asarray(data.mean(axis=0)))
    fano = gene_var/gene_mean
    top_genes = [i[0] for i in sorted(enumerate(gene_mean),\
        key=lambda x:x[1], reverse=True)][0:50] # index of top genes by mean
        
    A = (np.sqrt(gene_var)/gene_mean)[top_genes].min()
    B = 0.98
    expected_fano = (A**2)*gene_mean + (B**2)
        
    fstat = fano/expected_fano
    
    data_out = data[:,fstat>fthresh]
    
    # create boolean vector of genes to keep
    keep_genes = np.array([False]*gbm.matrix.shape[1])
    kg = np.array([False]*data.shape[1])
    kg[fstat>fthresh] = True
    keep_genes[cell_count>=cthresh] = kg
    
    # visualize Fano results
    plot_expected = np.array([[x,y] for (x,y) in sorted(zip(gene_mean,expected_fano))])
    
    if do_plot.lower()=='True'.lower():
        plt.figure()
        s=2   
        plt.scatter(gene_mean[fstat>fthresh],fano[fstat>fthresh],color='r',s=s)
        plt.scatter(gene_mean[fstat<=fthresh],fano[fstat<=fthresh],color='k',s=s)
        plt.plot(plot_expected[:,0],plot_expected[:,1])
        plt.title('Fano factor plot: %i high variance genes' %sum(fstat>fthresh))
        plt.xlabel('mu')
        plt.ylabel('Fano')
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.savefig('Fano.eps', format='eps')
        
    # plot number of cells with >count_thresh counts for all significant genes
        count_thresh = [1,2,5]
        Lct = len(count_thresh)
        ncells,ngenes = data_out.shape
        plot_col=2

        plt.figure()
        for i,ct in enumerate(count_thresh):
            counted_cells = np.sort((data_out>=ct).sum(axis=0))
            counted_cells = np.squeeze(np.asarray(counted_cells))[::-1] 
            proportion = counted_cells/float(ncells)
            plt.subplot(math.ceil(float(Lct)/plot_col), plot_col, i+1)
            plt.tight_layout()
            plt.plot(proportion,c='blue')
            plt.title('>= %i' %ct)
            plt.xlabel('genes')
            plt.ylabel('p')
        plt.savefig('sufficient_count_cells.eps',bbox_inches='tight',\
            format='eps',dpi=100)
        
    sys.stderr.write('done\n') 
    return keep_genes
# ------------------------------------------------------------------------------  
# force include genes
def force_include_genes(gbm,include_genes):
    force_keep_genes = np.array([False]*gbm.matrix.shape[1])
    for g in include_genes:
        try:
            index = np.where(gbm.ensembl==gbm.s2e[g])[0][0]
            force_keep_genes[index] = True
        except KeyError:
            print 'Gene %s does not exist'%g
        except IndexError:
            print 'Gene %s has zero count'%g  
    return force_keep_genes           
#-------------------------------------------------------------------------------
# Gene count normalizations

def normalize(data):
    sys.stderr.write('Normalizing... ')
    data_log = np.log2(data+1)
    sys.stderr.write('any nan: %s... '%np.any(np.isnan(data_log)))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
        data_norm = sp.stats.zscore(data_log,axis=0)
    data_norm = np.nan_to_num(data_norm)
    sys.stderr.write('done\n')
    return data_norm, data_log   
    
#-------------------------------------------------------------------------------
# Plot QC

def QCplot(data):
    sys.stderr.write('QC plots\n')
    """
    # UMI per cell
    UMI_per_cell = data.sum(axis=1)
    plt.figure()
    plt.hist(UMI_per_cell, bins=50)
    plt.ylabel('number of cells')
    plt.xlabel('UMI count')
    plt.savefig('UMI_per_cell.eps',bbox_inches='tight',\
        format='eps',dpi=100)
    """
      
#-------------------------------------------------------------------------------
# main   
def main(dir_, input_file, out_dir, size_factor_method, \
    thresholds, genome, outliers_file, include_genes, do_plot):
    
    os.chdir(out_dir)
    
    for i,t in enumerate(thresholds):
        if i==0:
            MT_cutoff = t
        if i==1:
            gpc_cutoff = t
        if i==2:
            UMIpc_cutoff = t    
        if i==3:
            erthresh = t
        if i==4:
            cthresh = t
        if i==5:    
            fthresh = t
            
    print '---- Loading gbm ----'       
    gbm_file = get_gbm_file(dir_, input_file, genome)
    print('File used: %s'%gbm_file)
    gbm = load_gbm(gbm_file,genome)
    print 'Matrix shape: %s '%(gbm.matrix.shape,)
    
    print '\n---- Removing genes with no counts ----'
    not_zero = gbm.matrix.getnnz(0)>0
    gbm[:,not_zero]
    print 'Matrix shape: %s'%(gbm.matrix.shape,)
        
    print '\n---- Removing bad cells ----'
    keep_cells = find_badcells(gbm, MT_cutoff, gpc_cutoff, UMIpc_cutoff, do_plot)
    gbm[keep_cells]
    print 'Matrix shape: %s'%(gbm.matrix.shape,) 
    
    print '\n---- Removing outliers ----'
    keep_cells = remove_outliers(gbm,outliers_file)
    gbm[keep_cells]
    print 'Matrix shape: %s'%(gbm.matrix.shape,)
    
    print '\n---- QC plots ----'
    QCplot(gbm.matrix)
    
    print '\n---- Adjusting for size factor and extracting signal ----'
    s = size_factor(gbm,size_factor_method, erthresh)
    keep_genes = get_signal(gbm,cthresh,fthresh,do_plot)
    force_keep_genes = force_include_genes(gbm,include_genes)
    GBM = copy.deepcopy(gbm) # non-reduced gbm
    gbm[:,keep_genes | force_keep_genes]
    print 'Matrix shape: %s'%(gbm.matrix.shape,)
    gbm.raw_counts = gbm.matrix
    GBM.raw_counts = GBM.matrix
    print 'Raw count matrix shape: %s'%(gbm.raw_counts.shape,)
    gbm.matrix = (gbm.matrix.T/s.T).T # convert to dense here
    GBM.matrix = (GBM.matrix.T/s.T).T
    
    print '\n---- Normalizing matrix ----'
    gbm.matrix, gbm.log = normalize(gbm.matrix)
    GBM.matrix, GBM.log = normalize(GBM.matrix)
    print 'Final check on reduced gbm components:'
    print 'Matrix shape: %s'%(gbm.matrix.shape,)
    print '%s'%type(gbm.matrix)
    print 'Barcodes shape: %s'%(gbm.barcodes.shape,)
    print '%s'%type(gbm.barcodes)
    print 'Gene list shape: %s'%(gbm.ensembl.shape,)
    print '%s'%type(gbm.ensembl)

    print 'Final check on full GBM components:'
    print 'Matrix shape: %s'%(GBM.matrix.shape,)
    print '%s'%type(GBM.matrix)
    print 'Barcodes shape: %s'%(GBM.barcodes.shape,)
    print '%s'%type(GBM.barcodes)
    print 'Gene list shape: %s'%(GBM.ensembl.shape,)
    print '%s'%type(GBM.ensembl)
    
    
    return gbm,GBM

if __name__ == '__main__':
    dir_ = os.getcwd()
    input_file = ''
    out_dir = dir_
    size_factor_method = 'total' # 'median' (DESeq) OR 'total' (TPM)
    thresholds = [0.12,2000,2000,0.01,10,1.5] # MT_cutoff, gpc_cutoff, UMIpc_cutoff erthresh, cthresh, fthresh
    genome = 'GRCh38'
    do_plot = 'True' # 'True' OR 'False' (case insensitive)
    markers = 'SOX2'
    outliers_file = ''
    include_genes = []
    
    ARG = sys.argv
    
    for i,arg in enumerate(ARG):
        if arg == '-f':
            input_file = ARG[i+1]
        if arg == '-d':
            dir_ = ARG[i+1]
        if arg == '-o':
            out_dir = ARG[i+1]
        if arg == '-s':
            size_factor_method = ARG[i+1]
        if arg == '-t':
            T = ARG[i+1].split(",")
            for i,el in enumerate(T):
                if el!='':
                    thresholds[i] = float(el)
        if arg == '-g':
            genome = ARG[i+1]
        if arg == '-p':
            do_plot = ARG[i+1]
        if arg == '-m':
            markers = ARG[i+1].split(",")
        if arg == '-z':
            outliers_file = ARG[i+1]
        if arg == '-i':
            include_genes = ARG[i+1].split(",")
        if arg == '-h':
            sys.stderr.write(Help)
            sys.exit(0)
            
            
    
    #data, data_count, s, ensembl2symb, symb2ensembl = \
    gbm,GBM = main(dir_, input_file, out_dir, size_factor_method, \
        thresholds, genome, outliers_file, include_genes, do_plot)
    
    # save data
    print '\n---- Saving data and making early visualization plots ----'
    do_save=True
    if do_save:
        """
        store = pd.HDFStore('gbm_adj.h5')
        store['data_norm'] = data
        store['data_count'] = data_count
        store['s'] = s
        store.close()
        """
        genes = np.array( [gbm.ensembl,np.array(['']*len(gbm.ensembl))],\
            dtype=np.object ).T
        for i in range(genes.shape[0]):
            genes[i,1] = gbm.e2s[genes[i,0]]

        GENES = np.array( [GBM.ensembl,np.array(['']*len(GBM.ensembl))],\
            dtype=np.object ).T
        for i in range(GENES.shape[0]):
            GENES[i,1] = GBM.e2s[GENES[i,0]]

        sio.savemat('gbm.mat',\
            {'matrix':gbm.matrix, 'logm':gbm.log, 'raw_counts':gbm.raw_counts, \
                'genes':genes, 'barcodes':gbm.barcodes, \
                'MATRIX':GBM.matrix, 'LOGM':GBM.log, 'RAW_COUNTS':GBM.raw_counts, \
                'GENES':GENES})
    
    # calculate pca
    if do_plot.lower()=='True'.lower():
        sys.stderr.write('Running pca... ')
        pca = PCA(n_components=2)
        pca.fit(gbm.matrix)
        Y = pca.transform(gbm.matrix)
        sys.stderr.write('done\n')
    
        # plot pca
        s=1
        plt.figure()
        plt.subplot(111)
        annotation = np.squeeze(np.asarray(GBM.raw_counts.sum(axis=1)))
        plt.scatter(Y[:,0],Y[:,1],c=annotation,edgecolors='face',\
            cmap='viridis',s=s)
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.colorbar()
        plt.savefig('pca.eps', format='eps')
        
        for m in markers:
            try:
                plt.figure()
                plt.subplot(111)
                annotation = np.squeeze(np.asarray(\
                    GBM.matrix[:,np.where(GBM.ensembl==GBM.s2e[m])[0][0]] ))
                plt.scatter(Y[:,0],Y[:,1],c=annotation,edgecolors='face',\
                    cmap='viridis',s=s)
                plt.xlabel('PC1')
                plt.ylabel('PC2')
                plt.colorbar()
                plt.savefig('%s.eps'%m, format='eps')
            except KeyError:
                print 'Marker %s not found in dataset'%m
            except IndexError:
                print 'Marker %s not found in dataset'%m
    
    
