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
from infra import *
################################################# (1) eucalidian distance between patients according tested genes expression profile #############################################

# (1) main
def find_expression_similarity_profile(gene_list_file_name, gene_expression_file_name, phenotype_file_name, gene_list_path=None, gene_expression_path=None, phenotype_path=None, source="GDC-TCGA",dataset="melanoma"):
    groups = load_expression_profile_by_labelling(gene_list_file_name, gene_expression_file_name, phenotype_file_name, tested_gene_path=gene_list_path, gene_expression_path=gene_expression_path, phenotype_path=phenotype_path, source=source, dataset=dataset)

    primary_labeled = groups[0]
    metastatic_labeled = groups[1]

    similarity_primary = []
    similarity_metastatic = []
    inter_similarity = []
    counter = 0
    for cur_1 in primary_labeled[1:]:
        counter += 1
        for cur_2 in primary_labeled[1:]:
            similarity_primary.append(np.linalg.norm(np.array(cur_1[1:], dtype=float)-np.array(cur_2[1:], dtype=float)))

    counter=0
    for cur_1 in metastatic_labeled[1:]:
        counter+=1
        for cur_2 in metastatic_labeled[1:]:
            similarity_metastatic.append(np.linalg.norm(np.array(cur_1[1:], dtype=float)-np.array(cur_2[1:], dtype=float)))

    counter = 0
    for cur_1 in primary_labeled[1:]:
        counter += 1
        for cur_2 in metastatic_labeled[1:]:
            inter_similarity.append(np.linalg.norm(np.array(cur_1[1:], dtype=float)-np.array(cur_2[1:], dtype=float)))

    print "similarity_primary: {}, {}".format(np.average(similarity_primary), np.var(similarity_primary))
    print "similarity_metastatic: {}, {}".format(np.average(similarity_metastatic), np.var(similarity_metastatic))
    print "inter_similarity: {}, {}".format(np.average(inter_similarity), np.var(inter_similarity))
    return np.average(similarity_primary), np.var(similarity_primary), np.average(similarity_metastatic), np.var(similarity_metastatic), np.average(inter_similarity), np.var(inter_similarity)


