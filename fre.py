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
import scipy.special
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from scipy.stats import hypergeom
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from svm import *
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

RANDOMIZED = "randomized"
REVERSED = "reversed"
NORMAL = "normal"

########################### FRE ###################################

def filter_gene_expressions(genelist_dataset, filter_ids_list, total_ids_list):
    included_genes_by_index = [i for i, gene_id in enumerate(total_ids_list) if gene_id in filter_ids_list]
    genelist_dataset = np.rot90(np.flip(genelist_dataset, 1), k=-1, axes=(1, 0))
    genelist_dataset = [cur for i, cur in enumerate(genelist_dataset) if i in included_genes_by_index]
    genelist_dataset = np.flip(np.rot90(genelist_dataset, k=1, axes=(1, 0)), 1)
    return genelist_dataset

def filter_gene_expressions_preprocessed(genelist_dataset, filter_ids_list, recursion_number_of_steps, recursion_step_size, start_index, total_ids_list):
    included_genes_by_index = [i for i, gene_id in enumerate(total_ids_list) if gene_id in filter_ids_list]
    included_genes_buckets = []
    included_matrices_prefixes = []
    for i in range(recursion_number_of_steps):
        included_genes_buckets.append([])
    genelist_dataset = np.rot90(np.flip(genelist_dataset, 1), k=-1, axes=(1, 0))
    for i, cur in enumerate(genelist_dataset):
        if i in included_genes_by_index:
            bucket_index = filter_ids_list.index(total_ids_list[i]) / recursion_step_size
            if bucket_index < len(included_genes_buckets):
                included_genes_buckets[bucket_index].append(cur)
    included_matrices_prefixes.append(included_genes_buckets[0])
    for cur in included_genes_buckets[1:]:
        included_matrices_prefixes.append(included_matrices_prefixes[-1] + cur)
    for i, cur in enumerate(included_matrices_prefixes):
        print i
        included_matrices_prefixes[i] = np.flip(np.rot90(cur, k=1, axes=(1, 0)), 1).tolist()
    return included_matrices_prefixes

def print_fre_results(test_scores, rounds, tested_gene_list_file_name, rank_method, permutation_type):
    print "about to print results for {} with rank method {} and permutation type {}".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')], rank_method, permutation_type)
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
def RFE(tested_gene_list_file_name, expression_profile_file_name, phenotype_file_name, rank_method=LOGISTIC_REGRESSION, gene_filter_file_name="protein_coding.txt", rounds=2, recursion_step_size=2, start_index = 0, recursion_number_of_steps=20, pval_preprocessing_file_name = None, permutation=NORMAL, groups=None, classification_method="svm_rbf_default", tuning_parameters={'C': [10], 'kernel': ['rbf']}):
    thismodule = sys.modules[__name__]
    clf = getattr(thismodule, classification_method)(tuning_parameters)
    print "about ot analyse: {}".format(tested_gene_list_file_name)
    # fetch gene expression by gene_id, divided by tumor type
    # test pval for significance differentiation between label values (primar vs metastatic)
    data, labels, groups, gene_ids = load_svm_data(tested_gene_list_file_name, expression_profile_file_name, phenotype_file_name,
                                                   gene_filter_file_name=gene_filter_file_name, groups = groups)
    if os.path.isfile(os.path.join(CACHE_DIR, pval_preprocessing_file_name)) and USE_CACHE:
        gene_pval_pair = load_sets(os.path.join(CACHE_DIR, pval_preprocessing_file_name))
        print "pval loaded from file"
    else:
        group_0_expression = groups[0]
        group_1_expression = groups[1]
        pvals = []
        gene_symbols = []
        for i in range(1, len(group_0_expression)):
            cur_pval = scipy.stats.ttest_ind([float(c) for c in group_0_expression[i][1:]],
                                             [float(c) for c in group_1_expression[i][1:]])[1]
            if not math.isnan(cur_pval):
                pvals.append(cur_pval)
                gene_symbols.append(group_0_expression[i][0])

        # sort gene_id-pval pairs by pval
        gene_pval_pair = zip(gene_symbols, pvals)
        gene_pval_pair.sort(key=lambda x: x[1], reverse=False)
        save_sets(gene_pval_pair,
                  os.path.join(CACHE_DIR, os.path.join(CACHE_DIR, pval_preprocessing_file_name)))
        print "pval saved to file"

    # calculate number of true hyphothesis after correction
    pvals = [cur[1] for cur in gene_pval_pair]
    fdr_results = fdrcorrection0(pvals, alpha=0.05, method='indep', is_sorted=True)
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    print "true hypothesis: {}/{}".format(true_counter, np.size(fdr_results[0]))

    gene_ids_ranked = [cur[0] for cur in gene_pval_pair]
    gene_ids_ranked = gene_ids_ranked[:true_counter]
    if permutation == RANDOMIZED:
        random.shuffle(gene_ids_ranked)
    elif permutation == REVERSED:
        gene_ids_ranked = list(reversed(gene_ids_ranked))  # random.shuffle(gene_ids_ranked)
    train_scores = []
    test_auPR_scores = []
    test_auROC_scores = []
    for j in range(recursion_number_of_steps):
        train_scores.append([])
        test_auPR_scores.append([])
        test_auROC_scores.append([])
    genelist_datasets = filter_gene_expressions_preprocessed(data, gene_ids_ranked, recursion_number_of_steps,
                                                             recursion_step_size, start_index, gene_ids)
    for i in range(rounds):
        genelist_datasets = np.rot90(genelist_datasets, k=1, axes=(1, 0))
        genelist_datasets, labels = randonize_patients(genelist_datasets, labels)
        genelist_datasets = np.rot90(genelist_datasets, k=-1, axes=(1, 0))
        for j in range(recursion_number_of_steps):
            # cur_dataset = filter_gene_expressions(genelist_dataset, gene_ids_ranked[:recursion_step_size*(j+1)], gene_ids)
            cur_dataset = genelist_datasets[j]
            data_train, data_test, labels_train, labels_test = divide_train_and_test_groups(cur_dataset, labels)
            test_auPR, test_auROC = apply_svm(clf, data_train, labels_train, data_test, labels_test, rank_method)
            test_auPR_scores[j].append(test_auPR)
            test_auROC_scores[j].append(test_auROC)
    print "#######################################"
    print "AUPR results:"
    avgs, vars = print_fre_results(test_auPR_scores, float(rounds), tested_gene_list_file_name, rank_method, permutation)
    print avgs
    print "AUROC results:"
    avgs, vars = print_fre_results(test_auROC_scores, float(rounds), tested_gene_list_file_name, rank_method, permutation)
    print avgs
