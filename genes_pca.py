import Bio.UniProt.GOA as GOA
# from orangecontrib.bio.go import Ontology
import wget
from utils.ensembl2gene_symbol import e2g_convertor
import time
import requests
import scipy.special
import matplotlib.pyplot as plt
from matplotlib import style
import matplotlib.ticker as ticker
style.use("ggplot")
import scipy
from lifelines.statistics import logrank_test
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
import openpyxl
from openpyxl import Workbook
from openpyxl.styles import Color, PatternFill, Font, Border, Side, Alignment
import sys
import os
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_ncbi_gene2go
# from goatools.associations import read_gaf
from lifelines import KaplanMeierFitter
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import rankdata
from genes_statistic import plot_genes_statistic
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
############################ () pca #############################


# () main
def genes_pca(tested_gene_list_file_name, total_gene_list_file_name, gene_expression_file_name, phenotype_file_name, survival_file_name, gene_filter_file_name=None, tested_gene_list_path=None, total_gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_file_path=None, var_th_index=None, is_unsupervised=True, start_k=2, end_k=2, meta_groups=None, filter_expression=None):
    data = load_integrated_ge_data(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, meta_groups=meta_groups, filter_expression=filter_expression)
    if data is None:
        print "insufficient data"
        return
    gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns, labels_assignment, survival_dataset = data
    plot_pca(gene_expression_top_var, labels_assignment, gene_expression_top_var_headers_columns, tested_gene_list_file_name, meta_groups)

def plot_pca(gene_expression_top_var, labels_assignment, gene_expression_top_var_headers_columns, tested_gene_list_file_name,meta_groups):
    actual_labels = list(range(1,len(meta_groups[0])+1))

    labels = [(cur["disease_code"]["value"][0], cur["_label"]) for i, cur in enumerate(meta_groups[0])]
    labels = [("unknown", 0)]+labels
    X=[]
    y=[]
    for cur in actual_labels:
        X.append(np.average(gene_expression_top_var[np.where(cur==labels_assignment[0]),:][0],axis=0))
        y.append(cur)
    X = np.array(X)
    y = np.array(y)
    pca = PCA(n_components=3)
    pca.fit_transform(X[:len(X)-1])

    fig = plt.figure(1, figsize=(20, 20))
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')
    # ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
    # y = np.choose(y, list(range(0, n_labels+1))).astype(np.float)
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
    for i, x in enumerate(X):
        name = labels[i+1]
        x2, y2, _ = proj3d.proj_transform(x[0], x[1], x[2], ax.get_proj())
        ax.annotate(name,
                    xy=(x2, y2), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))

    plt.savefig(os.path.join(constants.BASE_PROFILE, "output", "PCA_{}_{}.png").format(constants.CANCER_TYPE,time.time()))