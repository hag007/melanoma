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
from svm import *

########################### FRE ###################################

def filter_gene_expressions(genelist_dataset, filter_ids_list, total_ids_list):
    included_genes_by_index = [i for i, gene_id in enumerate(total_ids_list) if gene_id in filter_ids_list]
    genelist_dataset = np.rot90(np.flip(genelist_dataset, 1), k=-1, axes=(1, 0))
    genelist_dataset = [cur for i, cur in enumerate(genelist_dataset) if i in included_genes_by_index]
    genelist_dataset = np.flip(np.rot90(genelist_dataset, k=1, axes=(1, 0)), 1)
    return genelist_dataset

def filter_gene_expressions_preprocessed(genelist_dataset, filter_ids_list, recursion_number_of_steps, recursion_step_size, total_ids_list):
    included_genes_by_index = [i for i, gene_id in enumerate(total_ids_list) if gene_id in filter_ids_list]
    included_genes_buckets = []
    included_matrices_prefixes = []
    for i in range(recursion_number_of_steps):
        included_genes_buckets.append([])
    genelist_dataset = np.rot90(np.flip(genelist_dataset, 1), k=-1, axes=(1, 0))
    for i, cur in enumerate(genelist_dataset):
        if i in included_genes_by_index:
            bucket_index = filter_ids_list.index(total_ids_list[i])/recursion_step_size
            if bucket_index < len(included_genes_buckets):
                included_genes_buckets[bucket_index].append(cur)
    included_matrices_prefixes.append(included_genes_buckets[0])
    for cur in included_genes_buckets[1:]:
        included_matrices_prefixes.append(included_matrices_prefixes[-1] + cur)
    for i, cur in enumerate(included_matrices_prefixes):
        included_matrices_prefixes[i] = np.flip(np.rot90(cur, k=1, axes=(1, 0)), 1).tolist()
    return included_matrices_prefixes

def print_fre_results(test_scores, rounds):
    avgs = []
    vars = []
    for i, cur_i in enumerate(test_scores):
        avg = sum(test_scores[i]) / rounds
        var = np.var(test_scores[i])
        print "group {}, train_accuracy_avg: {}, train_accuracy_var: {}".format(i, avg, var)
        avgs.append(avg)
        vars.append(var)
    return avgs, vars

# (4) main
def RFE(tested_gene_list_file_name, expression_profile_file_name, phenotype_file_name, gene_filter_file_name="protein_coding.txt", rounds=2, recursion_step_size=2, recursion_number_of_steps=20, pval_preprocessing_file_name = None):

    print "about ot analyse: {}".format(tested_gene_list_file_name)
    # fetch gene expression by gene_id, divided by tumor type
    data, labels, primary_expression, metastatic_expression, gene_ids = load_svm_data(tested_gene_list_file_name, expression_profile_file_name, phenotype_file_name,
                                 gene_filter_file_name)
    # test pval for significance differentiation between label values (primar vs metastatic)
    if os.path.isfile(os.path.join(BASE_OUTPUT_DIR, pval_preprocessing_file_name)) and IS_FORCE == False:
        gene_pval_pair = load_sets(os.path.join(BASE_OUTPUT_DIR, pval_preprocessing_file_name))
        print "pval loaded from file"
    else:
        primary_expression = np.rot90(np.flip(primary_expression, 1), k=-1, axes=(1, 0))
        metastatic_expression = np.rot90(np.flip(metastatic_expression, 1), k=-1, axes=(1, 0))
        pvals = []
        gene_symbols = []
        for i in range(1, len(primary_expression)):
            cur_pval = scipy.stats.ttest_ind([float(c) for c in primary_expression[i][1:]],
                                             [float(c) for c in metastatic_expression[i][1:]])[1]
            if not math.isnan(cur_pval):
                pvals.append(cur_pval)
                gene_symbols.append(primary_expression[i][0])

        # sort gene_id-pval pairs by pval
        gene_pval_pair = zip(gene_symbols, pvals)
        gene_pval_pair.sort(key=lambda x: x[1], reverse=False)
        save_sets(gene_pval_pair,
                  os.path.join(BASE_OUTPUT_DIR, os.path.join(BASE_OUTPUT_DIR, pval_preprocessing_file_name)))
        print "pval saved to file"

    gene_ids_ranked = [cur[0] for cur in gene_pval_pair]

    train_scores = []
    test_scores = []
    for j in range(recursion_number_of_steps):
        train_scores.append([])
        test_scores.append([])
    genelist_datasets = filter_gene_expressions_preprocessed(data, gene_ids_ranked, recursion_number_of_steps,
                                                             recursion_step_size, gene_ids)
    for i in range(rounds):
        genelist_datasets = np.rot90(genelist_datasets, k=1, axes=(1, 0))
        genelist_datasets, labels = randonize_patients(genelist_datasets, labels)
        genelist_datasets = np.rot90(genelist_datasets, k=-1, axes=(1, 0))
        for j in range(recursion_number_of_steps):
            # cur_dataset = filter_gene_expressions(genelist_dataset, gene_ids_ranked[:recursion_step_size*(j+1)], gene_ids)
            cur_dataset = genelist_datasets[j]
            data_train, data_test, labels_train, labels_test = divide_train_and_test_groups(cur_dataset, labels)
            tuned_parameters = {'C': [10], 'kernel': ['rbf']}
            test_accuracy = apply_svm(tuned_parameters, data_train, labels_train, data_test, labels_test)
            test_scores[j].append(test_accuracy)
    print "#######################################"
    avgs, vars = print_fre_results(test_scores, float(rounds))
    print avgs
