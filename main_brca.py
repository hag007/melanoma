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
from gene_sets_correlation import *
#from svm import *
import fre
from fre import *
from feature_selection import *
from utils.param_builder import *
from genes_clustering import *


##################### BRCA ########################
# update_dirs(DATASET_DIR="TCGA\\breast\\")

## prediction
# prediction_by_gene_expression(["mito.txt", "apoptosis_genes.txt", "lipid_synthase_genes.txt", "oxidative.txt", "lysosome_genes.txt", "oxidation_genes.txt", "top420s.txt"], "TCGA-BRCA.htseq_fpkm.tsv", "BRCA_clinicalMatrix", label="breast_carcinoma_estrogen_receptor_status", label_values=["Positive","Negative"], rank_method=LOGISTIC_REGRESSION, gene_filter_file_name="protein_coding.txt", rounds=30)

## randomized by size
# RFE("protein_coding.txt", "TCGA-BRCA.htseq_fpkm.tsv", "BRCA_clinicalMatrix", label="PAM50Call_RNAseq", label_values=["LumA", "LumB"], rank_method = LOGISTIC_REGRESSION, rounds=30, recursion_step_size=13, recursion_number_of_steps=1, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=RANDOMIZED)
# RFE("protein_coding.txt", "TCGA-BRCA.htseq_fpkm.tsv", "BRCA_clinicalMatrix", label="PAM50Call_RNAseq", label_values=["LumA", "LumB"], rank_method = LOGISTIC_REGRESSION, rounds=30, recursion_step_size=120, recursion_number_of_steps=1, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=RANDOMIZED)
# RFE("protein_coding.txt", "TCGA-BRCA.htseq_fpkm.tsv", "BRCA_clinicalMatrix", label="PAM50Call_RNAseq", label_values=["LumA", "LumB"], rank_method = LOGISTIC_REGRESSION, rounds=30, recursion_step_size=420, recursion_number_of_steps=1, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=RANDOMIZED)
# RFE("protein_coding.txt", "TCGA-BRCA.htseq_fpkm.tsv", "BRCA_clinicalMatrix", label="PAM50Call_RNAseq", label_values=["LumA", "LumB"], rank_method = LOGISTIC_REGRESSION, rounds=30, recursion_step_size=750, recursion_number_of_steps=1, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=RANDOMIZED)

# RFE("protein_coding.txt", "TCGA-BRCA.htseq_fpkm.tsv", "BRCA_clinicalMatrix", label="metastatic_breast_carcinoma_estrogen_receptor_status", label_values=["Positive","Negative"], rank_method = LOGISTIC_REGRESSION, rounds=30, recursion_step_size=13, recursion_number_of_steps=1, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=RANDOMIZED)
# RFE("protein_coding.txt", "TCGA-BRCA.htseq_fpkm.tsv", "BRCA_clinicalMatrix", label="metastatic_breast_carcinoma_estrogen_receptor_status", label_values=["Positive","Negative"], rank_method = LOGISTIC_REGRESSION, rounds=30, recursion_step_size=120, recursion_number_of_steps=1, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=RANDOMIZED)
# RFE("protein_coding.txt", "TCGA-BRCA.htseq_fpkm.tsv", "BRCA_clinicalMatrix", label="metastatic_breast_carcinoma_estrogen_receptor_status", label_values=["Positive","Negative"], rank_method = LOGISTIC_REGRESSION, rounds=30, recursion_step_size=420, recursion_number_of_steps=1, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=RANDOMIZED)
# RFE("protein_coding.txt", "TCGA-BRCA.htseq_fpkm.tsv", "BRCA_clinicalMatrix", label="metastatic_breast_carcinoma_estrogen_receptor_status", label_values=["Positive","Negative"], rank_method = LOGISTIC_REGRESSION, rounds=30, recursion_step_size=750, recursion_number_of_steps=1, pval_preprocessing_file_name = "pvals_protein_coding.txt", permutation=RANDOMIZED)

