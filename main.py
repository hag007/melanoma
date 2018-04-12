import re
import os
import math
from numpy import *
from decimal import Decimal, ROUND_HALF_EVEN
from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction
import numpy.random
from sklearn.datasets import fetch_mldata
import sklearn.preprocessing
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import scipy.special
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
from sklearn import svm
from sklearn import svm
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.model_selection import PredefinedSplit
from sklearn.metrics import accuracy_score
import scipy
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import time
from matplotlib.ticker import FormatStrFormatter
import math
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from clusters_similarity import *
from significant_enrichment import *
from svm import *
from fre import *

# (1) similarity caller:
# find_expression_similarity_profile(gene_list_file_name="oxidative.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", phenotpe_file_name="TCGA-SKCM.GDC_phenotype.tsv")
# (2) call:
find_expression_significance("mito.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_8490_13.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=8490, B = 13)
# find_expression_significance("lysosome_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_19976_754.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=19976, B = 754)
# find_expression_significance("apoptosis_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_19976_119.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=19976, B = 119)
# find_expression_significance("lipid_synthase_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_19976_129.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=19976, B = 129)
# find_expression_significance("oxidation_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_19976_1066.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=19976, B = 1066)
# find_expression_significance("top420s.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_19976_420.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=19976, B = 420)
# find_expression_significance("oxidative.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_19976_420.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=19976, B = 420)
#------ for FDR results only:
# find_expression_significance("mito.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_8490_13.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=8490, B = 13)
# find_expression_significance("lysosome_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_8490_329.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=8490, B = 329)
# find_expression_significance("apoptosis_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_8490_59.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=8490, B = 59)
# find_expression_significance("lipid_synthase_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_8490_53.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=8490, B = 53)
# find_expression_significance("oxidation_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_8490_421.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=8490, B = 421)
# find_expression_significance("top420s.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_8490_421.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=8490, B = 420)
# find_expression_significance("oxidative.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_8490_194.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=8490, B = 194)
# (3) call
# predict_tumor_type_by_tested_gene_expression(["oxidative.txt", "oxidation_genes.txt", "mito.txt", "apoptosis_genes.txt", "lysosome_genes.txt", "lipid_synthase_genes.txt", "top420s.txt"], "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rounds=30)

# (4)
# RFE("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rounds=2, recursion_step_size=5, recursion_number_of_steps=100, pval_preprocessing_file_name = "pvals_protein_coding.txt")


# primaries_avg = []
# primaries_var = []
# metstatics_avg = []
# metstatics_var = []
# inters_avg = []
# inters_var = []
#
# for i in range(100):
#     print "interation {}".format(i)
#     primary_similarity_avg, primary_similarity_var, metastatic_similarity_avg, metastatic_similarity_var, inter_similarity_avg, inter_similarity_var = find_expression_similarity_profile()
#     print "primary:{}, {}; metastatic:{}, {}; inter: {}, {}".format(primary_similarity_avg, primary_similarity_var, metastatic_similarity_avg, metastatic_similarity_var, inter_similarity_avg, inter_similarity_var)
#
#
#     primaries_avg.append(primary_similarity_avg)
#     primaries_var.append(primary_similarity_var)
#     metstatics_avg.append(metastatic_similarity_avg)
#     metstatics_var.append(metastatic_similarity_var)
#     inters_avg.append(inter_similarity_avg)
#     inters_var.append(inter_similarity_var)
#
# print primaries_avg
# print primaries_var
# print metstatics_avg
# print metstatics_var
# print inters_avg
# print inters_var


