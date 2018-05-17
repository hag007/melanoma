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
# import matplotlib as mpl
#mpl.use('Agg')
import scipy.special
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
# import matplotlib.pyplot as plt
# from matplotlib import style
# style.use("ggplot")
from sklearn import svm
from sklearn import svm
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.model_selection import PredefinedSplit
from sklearn.metrics import accuracy_score
import scipy
from scipy.stats import hypergeom
# from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import time
# from matplotlib.ticker import FormatStrFormatter
import math
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
#from clusters_similarity import *
from significant_enrichment import *
#from svm import *
from fre import *

# (1) similarity caller:
# find_expression_similarity_profile(gene_list_file_name="oxidative.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", phenotpe_file_name="TCGA-SKCM.GDC_phenotype.tsv")
# (2) call:
# find_expression_significance("mito.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_19976_13.npy", pval_preprocessing_file_name = "pvals_protein_coding.txt", N=19976, B = 13)
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

#------ mirna ------------:
# find_expression_significance("mir-MT-ATP6-244.txt", "mir_total.txt", "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_1610_244.npy", pval_preprocessing_file_name = "pvals_mirna.txt", N=1610, B = 244)
# find_expression_significance("mir-MT-ATP8-161.txt", "mir_total.txt", "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_1610_161.npy", pval_preprocessing_file_name = "pvals_mirna.txt", N=1610, B = 161)
# find_expression_significance("mir-MT-ND4L-427.txt", "mir_total.txt", "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_1610_427.npy", pval_preprocessing_file_name = "pvals_mirna.txt", N=1610, B = 427)
# find_expression_significance("mir-MT-ND5-161.txt", "mir_total.txt", "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_1610_161.npy", pval_preprocessing_file_name = "pvals_mirna.txt", N=1610, B = 161)
# find_expression_significance("mir-mt-total.txt", "mir_total.txt", "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_1610_707.npy", pval_preprocessing_file_name = "pvals_mirna.txt", N=1610, B = 707)
# find_expression_significance("mir-shc1-493.txt", "mir_total.txt", "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", hgt_preprocessing_file_name = "HGTs_out_1610_493.npy", pval_preprocessing_file_name = "pvals_mirna.txt", N=1610, B = 493)

# mHGT_list = ['mir-tami_metabolism-ATP5G.txt', 'mir-tami_metabolism-SOD2.txt', 'mir-tami_metabolism-ATP5D.txt', 'mir-tami_metabolism-CATALASE.txt', 'mir-tami_metabolism-PPARGC1A.txt', 'mir-tami_metabolism-SLC25A4.txt', 'mir-tami_metabolism-PPRC1.txt', 'mir-tami_metabolism-NDUFA8.txt', 'mir-tami_metabolism-PCG1A.txt', 'mir-tami_metabolism-SOD1.txt', 'mir-tami_metabolism-ATP5B.txt', 'mir-tami_mito-SLC2A4.txt', 'mir-tami_mito-Hk2.txt', 'mir-tami_mito-PKM.txt', 'mir-tami_mito-Pdk1.txt', 'mir-tami_mito-LDHA.txt', 'mir-tami_extra-MCT4.txt', 'mir-tami_extra-MCM10.txt', 'mir-tami_extra-POLD3.txt', 'mir-tami_extra-SLC2A1.txt', 'mir-tami_extra-IGFBP3.txt', 'mir-tami_extra-hk1.txt']
# pvals = []
# for cur in mHGT_list:
#     pvals.append(find_expression_significance(cur, "mir_total.txt", "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", pval_preprocessing_file_name = "pvals_mirna.txt", N=1610))
# print mHGT_list
# print pvals


# (3) call
# predict_tumor_type_by_tested_gene_expression(["oxidative.txt", "oxidation_genes.txt", "mito.txt", "apoptosis_genes.txt", "lysosome_genes.txt", "lipid_synthase_genes.txt", "top420s.txt"], "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rank_method=DISTANCE, gene_filter_file_name="protein_coding.txt", rounds=100)
# predict_tumor_type_by_tested_gene_expression(["oxidative.txt", "oxidation_genes.txt", "mito.txt", "apoptosis_genes.txt", "lysosome_genes.txt", "lipid_synthase_genes.txt", "top420s.txt"], "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rank_method=LOGISTIC_REGRESSION, gene_filter_file_name="protein_coding.txt", rounds=100)
# predict_tumor_type_by_tested_gene_expression(["test.txt"], "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", gene_filter_file_name="protein_coding.txt", rounds=10)
predict_tumor_type_by_tested_gene_expression(["mito.txt"], "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rank_method=LOGISTIC_REGRESSION, gene_filter_file_name="protein_coding.txt", rounds=30)
predict_tumor_type_by_tested_gene_expression(["mito.txt"], "TCGA-SKCM.htseq_fpkm.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rank_method=LOGISTIC_REGRESSION, gene_filter_file_name="protein_coding.txt", rounds=30)
predict_tumor_type_by_tested_gene_expression(["mito.txt"], "TCGA-SKCM.htseq_fpkm-uq.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rank_method=LOGISTIC_REGRESSION, gene_filter_file_name="protein_coding.txt", rounds=30)
#------ mirna ------------:
# predict_tumor_type_by_tested_gene_expression(["mir-MT-ATP6-244.txt", "mir-MT-ATP8-161.txt", "mir-MT-ND4L-427.txt", "mir-MT-ND5-161.txt"], "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rank_method=DISTANCE, rounds=100)
# predict_tumor_type_by_tested_gene_expression(['mir-tami_metabolism-ATP5G.txt', 'mir-tami_metabolism-SOD2.txt', 'mir-tami_metabolism-ATP5D.txt', 'mir-tami_metabolism-CATALASE.txt', 'mir-tami_metabolism-PPARGC1A.txt', 'mir-tami_metabolism-SLC25A4.txt', 'mir-tami_metabolism-PPRC1.txt', 'mir-tami_metabolism-NDUFA8.txt', 'mir-tami_metabolism-PCG1A.txt', 'mir-tami_metabolism-SOD1.txt', 'mir-tami_metabolism-ATP5B.txt', 'mir-tami_mito-SLC2A4.txt', 'mir-tami_mito-Hk2.txt', 'mir-tami_mito-PKM.txt', 'mir-tami_mito-Pdk1.txt', 'mir-tami_mito-LDHA.txt', 'mir-tami_extra-MCT4.txt', 'mir-tami_extra-MCM10.txt', 'mir-tami_extra-POLD3.txt', 'mir-tami_extra-SLC2A1.txt', 'mir-tami_extra-IGFBP3.txt', 'mir-tami_extra-hk1.txt'], "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rank_method=DISTANCE, rounds=30)
# predict_tumor_type_by_tested_gene_expression(["mir-MT-ATP6-244.txt", "mir-MT-ATP8-161.txt", "mir-MT-ND4L-427.txt", "mir-MT-ND5-161.txt"], "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rank_method=LOGISTIC_REGRESSION, rounds=100)
# predict_tumor_type_by_tested_gene_expression(["mir_total.txt"], "TCGA-SKCM.mirna.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rank_method=LOGISTIC_REGRESSION, rounds=30)

# (4)
# RFE("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv",rank_method = DISTANCE, rounds=5, recursion_step_size=100, recursion_number_of_steps=80, start_index = 0, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=NORMAL)
# RFE("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv",rank_method = LOGISTIC_REGRESSION, rounds=5, recursion_step_size=50, recursion_number_of_steps=200, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=NORMAL)
# RFE("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv",rank_method = DISTANCE, rounds=5, recursion_step_size=50, recursion_number_of_steps=200, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=RANDOMIZED)
# RFE("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv",rank_method = LOGISTIC_REGRESSION, rounds=5, recursion_step_size=50, recursion_number_of_steps=200, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=RANDOMIZED)
# RFE("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv",rank_method = DISTANCE, rounds=2, recursion_step_size=200, recursion_number_of_steps=50, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=REVERSED)
# RFE("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv",rank_method = LOGISTIC_REGRESSION, rounds=2, recursion_step_size=200, recursion_number_of_steps=50, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=REVERSED)



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


