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


############################################(3) SVM to see whether we can predict tumor type by tested genes' profile expression########################################

def load_svm_data(tested_gene_file_name, expression_profile_file_name, phenotype_file_name, gene_filter_file_name=None):
    data = []
    labels = []
    primary_labeled, metastatic_labeled = load_expression_profile_by_gene_and_tumor_type(tested_gene_file_name, expression_profile_file_name, phenotype_file_name, gene_filter_file_name)

    for cur in primary_labeled[1:]:
        data.append(cur[1:].astype(np.float))
        labels.append(0)

    for cur in metastatic_labeled[1:]:
        data.append(cur[1:].astype(np.float))
        labels.append(1)

    return data, labels, primary_labeled, metastatic_labeled, primary_labeled[0]

def randonize_patients(data, labels):
    total = np.c_[data, labels]
    random.shuffle(total)
    # results = zip(*total)
    return total[:,:-1], total[:,-1]



def divide_train_and_test_groups(data, labels):
    data_train = data[:(3 * len(data)) / 4]
    data_test = data[(3 * len(data)) / 4:]
    labels_train = labels[:(3 * len(data)) / 4]
    labels_test = labels[(3 * len(data)) / 4:]
    return data_train, data_test, labels_train, labels_test

def apply_svm(tuned_parameters, data_train, labels_train, data_test, labels_test):
    gs_train = GridSearchCV(svm.SVC(probability=True), param_grid=tuned_parameters, return_train_score=True)
    data_train = [[cur2 for cur2 in cur1] for cur1 in data_train]
    labels_train = [i for i in labels_train]
    data_test = [[cur2 for cur2 in cur1] for cur1 in data_test]
    labels_test = [i for i in labels_test]
    gs_train.fit(data_train, labels_train)
    # predicted_results = gs_train.predict(data_test)
    probabilities = gs_train.decision_function(data_test)
    # probabilities = gs_train.predict_proba(data_test)
    # probabilities = probabilities[:,1]
    #zipped = zip(probabilities, data_train, labels_train)
    #zipped_sorted = sorted(zipped, key=lambda x: x[0])
    #patients_rank_sorted = [x[0] for x in zipped_sorted]
    #patients_expression_sorted = [x[1] for x in zipped_sorted]
    #patients_labels_sorted = [x[2] for x in zipped_sorted]
    # precision, recall, _ = precision_recall_curve(labels_test, probabilities)
    average_precision = average_precision_score(labels_test, probabilities)
    print "average_precision: {}".format(average_precision)
    ###
    #
    # plt.step(recall, precision, color='b', alpha=0.2,
    #          where='post')
    # plt.fill_between(recall, precision, step='post', alpha=0.2,
    #                  color='b')
    #
    # plt.xlabel('Recall')
    # plt.ylabel('Precision')
    # plt.ylim([0.0, 1.05])
    # plt.xlim([0.0, 1.0])
    # plt.title('2-class Precision-Recall curve: AP={0:0.2f}'.format(average_precision))
    # plt.savefig(os.path.join(BASE_OUTPUT_DIR,"test.png"))

    ###
    # ta = (len(labels_test) * 1.0 - sum([abs(p - r) for p, r in zip(predicted_results, labels_test)])) / (len(labels_test) * 1.0)
    # return gs_train.best_score_, ta
    return average_precision

def print_svm_results(train_scores, train_alt_scores, test_scores, test_alt_scores, rounds):
    print "train_accuracy_avg: {}".format(sum(train_scores)/ rounds)
    print "train_alt_accuracy_avg: {}".format(sum(train_alt_scores)/ rounds)
    print "train_accuracy_diff_avg: {}".format((sum(train_scores)-sum(train_alt_scores))/rounds)
    print "test_accuracy_avg: {}".format(sum(test_scores)/ rounds)
    print "test_alt_accuracy_avg: {}".format(sum(test_alt_scores)/ rounds)
    print "test_accuracy_diff_avg: {}".format((sum(test_scores)-sum(test_alt_scores))/rounds)
    print "train p val: {}".format(scipy.stats.ttest_ind(train_scores, train_alt_scores)[1])
    print "test_p_val: {}".format(scipy.stats.ttest_ind(test_scores, test_alt_scores)[1])

    # print "{}".format(train_score_sum / rounds)
    # print "{}".format(train_alt_score_sum / rounds)
    # print "{}".format(scipy.stats.ttest_ind(train_scores, train_alt_scores)[1])
    # print "{}".format(test_score_sum / rounds)
    # print "{}".format(test_alt_score_sum / rounds)
    # print "{}".format(scipy.stats.ttest_ind(test_scores, train_alt_scores)[1])

def print_svm_results(test_scores, rounds):
    avgs = []
    for i, cur_i in enumerate(test_scores):
        avg = sum(test_scores[i]) / rounds
        print "group {}, train_accuracy_avg: {}, train_accuracy_var: {}".format(i, avg, np.var(test_scores[i]))
        avgs.append(avg)

    results_summarized = []
    for i, cur_i in enumerate(avgs):
        results_summarized.append([])

    for i, cur_i in enumerate(avgs):
        for j, cur_j in enumerate(avgs):
            if j>i:
                pval = scipy.stats.ttest_ind(test_scores[i], test_scores[j])[1]
                print "test_p_val for {}, {}: {}".format(i,j,pval)
                results_summarized[i].append(pval)
            elif j==i:
                results_summarized[i].append(1)
            else:
                results_summarized[i].append(-1)

    for i, cur_i in enumerate(avgs):
        for j, cur_j in enumerate(avgs):
            if results_summarized[i][j] == -1:
                results_summarized[i][j] = results_summarized[j][i]

    for i, cur_i in enumerate(avgs):
        print ""
        for j, cur_j in enumerate(avgs):
            print "{}\t".format(results_summarized[i][j]),

    # print "{}".format(train_score_sum / rounds)
    # print "{}".format(train_alt_score_sum / rounds)
    # print "{}".format(scipy.stats.ttest_ind(train_scores, train_alt_scores)[1])
    # print "{}".format(test_score_sum / rounds)
    # print "{}".format(test_alt_score_sum / rounds)
    # print "{}".format(scipy.stats.ttest_ind(test_scores, train_alt_scores)[1])


# (3) main
def predict_tumor_type_by_tested_gene_expression(tested_gene_file_names, expression_profile_file_name, phenotype_file_name, gene_filter_file_name="protein_coding.txt", rounds=2):

    genelist_datasets = []
    for tested_gene_file_name in tested_gene_file_names:
        data, labels , _1, _2, _3 = load_svm_data(tested_gene_file_name, expression_profile_file_name, phenotype_file_name,
                                     gene_filter_file_name)
        genelist_datasets.append(data)

    train_scores = []
    test_scores = []
    for j in range(len(genelist_datasets)):
        train_scores.append([])
        test_scores.append([])
    for i in range(rounds):
        genelist_datasets = np.rot90(genelist_datasets, k=-1, axes=(1, 0))
        genelist_datasets, labels = randonize_patients(genelist_datasets, labels)
        genelist_datasets = np.rot90(genelist_datasets, k=1, axes=(1, 0))
        for j, cur_dataset in enumerate(genelist_datasets):
            data_train, data_test, labels_train, labels_test = divide_train_and_test_groups(cur_dataset, labels)
            tuned_parameters = {'C': [10], 'kernel': ['rbf']}
            test_accuracy = apply_svm(tuned_parameters, data_train, labels_train, data_test, labels_test)
            test_scores[j].append(test_accuracy)
    print "#######################################"
    print_svm_results(test_scores, float(rounds))
