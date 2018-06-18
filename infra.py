import time
import re
import os
import math
import random
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
import matplotlib.pyplot as plt
from utils.stopwatch import Stopwatch
from matplotlib import style
style.use("ggplot")
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
import constants
############################################ infra ########################################

def load_gene_list(gene_list_file_name, gene_list_path=None): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.LIST_DIR,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines

# return gene expression table filtered according an external list with proper orientation (0 degree angle according genes, 90 degree angle according patients)
def load_gene_expression_profile(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None ,by_gene=False):
    stopwatch = Stopwatch()
    stopwatch.start()
    gene_list = load_gene_list(gene_list_file_name=gene_list_file_name, gene_list_path=gene_list_path)
    gene_list = [l.split(".")[0] for i, l in enumerate(gene_list)]
    print stopwatch.stop("done loading gene list")
    # random.shuffle(gene_list)
    # gene_list = gene_list[:400]
    if gene_filter_file_name:
        stopwatch.start()
        filter_gene_list = load_gene_list(gene_list_file_name=gene_filter_file_name, gene_list_path=gene_filter_path)
        gene_list = [cur for cur in gene_list if cur in filter_gene_list]
        print stopwatch.stop("done filter gene list")

    if gene_expression_path == None:
        gene_expression_path = os.path.join(constants.TCGA_DATA_DIR, gene_expression_file_name)
        stopwatch.start()
    f = open(gene_expression_path,'r')
    expression_profiles_filtered = [l.strip().split() for i, l in enumerate(f) if i==0 or l[:l.strip().find('\t')].split(".")[0] in gene_list]
    # or l.strip()[0:l.strip().find('\t')] in gene_list or l.strip()[0:l.strip().find('\t')].split(".")[0] in gene_list
    f.close()
    print stopwatch.stop("done filter gene expression")
    if not by_gene:
        stopwatch.start()
        expression_profiles_filtered = np.flip(np.rot90(expression_profiles_filtered, k=1, axes=(1,0)),1)
        print stopwatch.stop("done rotate gene expression")

    return expression_profiles_filtered


def load_gene_expression_profile_by_genes(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None):
    return load_gene_expression_profile(gene_list_file_name=gene_list_file_name, gene_expression_file_name=gene_expression_file_name, gene_filter_file_name=gene_filter_file_name, gene_list_path=gene_list_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path, by_gene=True)

def load_gene_expression_profile_by_patients(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None):
    return load_gene_expression_profile(gene_list_file_name=gene_list_file_name, gene_expression_file_name=gene_expression_file_name, gene_filter_file_name=gene_filter_file_name, gene_list_path=gene_list_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path,by_gene=False)


def load_phenotype_data(phenotype_file_name, phenotype_list_path=None, source="GDC-TCGA",dataset="melanoma"):
    if not phenotype_list_path:
        phenotype_list_path = os.path.join(constants.TCGA_DATA_DIR,phenotype_file_name)
    f = open(phenotype_list_path, 'r')
    phenotype_profile = [l.strip().split('\t') for l in f]
    f.close()
    return phenotype_profile

def load_survival_data(survival_file_name, survival_list_path=None, source="GDC-TCGA",dataset="melanoma"):
    if not survival_list_path:
        phenotype_list_path = os.path.join(constants.TCGA_DATA_DIR,survival_file_name)
    f = open(phenotype_list_path, 'r')
    phenotype_profile = [l.strip().split('\t') for l in f]
    f.close()
    return phenotype_profile


def divided_patient_ids_by_label(phenotype_list_file_name, phenotype_list_path=None, labels=None, label_values=None, groups=None):
    if not groups and not labels:
        groups = [{"sample_type.samples" :{"type": "string", "value": ["Primary Tumor"]}},
                  {"sample_type.samples": {"type": "string", "value": ["Metastatic"]}}]
    elif not groups and any(labels):
        divided_patient_ids_by_label_old(phenotype_list_file_name, phenotype_list_path, labels, label_values)
        return
    phenotype_data_formatted = load_phenotype_data(phenotype_list_file_name, phenotype_list_path)
    headers = phenotype_data_formatted[0]
    phenotype_profiles = phenotype_data_formatted[1:]
    patients_by_labeling = []
    for i in range(len(groups)):
        patients_by_labeling.append([])
    label_indices = []
    duplicates_counter = 0
    try:
        for i, pp in enumerate(phenotype_profiles):
            for j, cur_group in enumerate(groups):
                dup=0
                is_hold_constraits = True
                for k,v in cur_group.iteritems():
                    if len(pp) <= headers.index(k): continue
                    if v['type'] == "string":
                        if not any([pp[headers.index(k)] == cur for cur in v["value"]]):
                            is_hold_constraits = False
                if is_hold_constraits:
                    patients_by_labeling[j].append(pp[0])
                    dup+=1
            if dup > 1:
                duplicates_counter+=1
        print "number of duplicated patients: {}".format(duplicates_counter)
    except ValueError:
        print ValueError
        pass
    return (patients_by_labeling)


def divided_patient_ids_by_label_old(phenotype_list_file_name, phenotype_list_path=None, labels=None, label_values=None,
                                 group_0=None, group_1=None):
    if not labels:
        labels = [constants.LABEL_ID]
    if type(labels) is not list:
        labels = [labels]
    if not label_values:
        label_values = [[constants.PRIMARY_TUMOR], [constants.METASTATIC]]
    if type(label_values[0]) is not list:
        label_values = [[label_values[0]], [label_values[1]]]
    phenotype_data_formatted = load_phenotype_data(phenotype_list_file_name, phenotype_list_path)
    headers = phenotype_data_formatted[0]
    phenotype_profiles = phenotype_data_formatted[1:]
    label_0_patients = []
    label_1_patients = []
    label_indices = []
    for label in labels:
        label_indices.append([i for i, v in enumerate(headers) if v == label][0])
    for pp in phenotype_profiles:
        if all([pp[cur_idx] == label_values[0][i] for i, cur_idx in enumerate(label_indices)]):
            label_0_patients.append(pp[0])
        elif all([pp[cur_idx] == label_values[1][i] for i, cur_idx in enumerate(label_indices)]):
            label_1_patients.append(pp[0])
    return (label_0_patients, label_1_patients)

# def transform_to_map(tbl):
#     map = {}
#     for cur in tbl[1:]:
#         map[cur[0]] = cur[1:]
#     return map

# load expression profile filtered by an external genes list and divided according tumor type label
def load_expression_profile_by_labelling(gene_list_file_name, gene_expression_file_name, phenotype_file_name, label=None, label_values=None, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_path=None, groups=None):

    expression_profiles_formatted = load_gene_expression_profile_by_patients(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=gene_filter_file_name, gene_list_path=gene_list_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path)
    patients_by_labeling = divided_patient_ids_by_label(phenotype_file_name, phenotype_path, label, label_values, groups)

    expression_profiles_by_labeling = []
    for i in patients_by_labeling:
        expression_profiles_by_labeling.append([])
    logger.info("about to split expression by primary tumor and metastatic")

    logger.info("expression_profile size: {},{}".format(*np.shape(expression_profiles_formatted)))
    for i,cur in enumerate(expression_profiles_formatted):
        if i==0: # that is, vertical headers
            for cur_group in expression_profiles_by_labeling:
                cur_group.append(cur)
        cur_id = expression_profiles_formatted[i][0]
        if not str.isdigit(cur_id[-1]) and constants.PHENOTYPE_FORMAT == "TCGA":
            cur_id = cur_id[:-1]
        patient_found = False
        for k, cur_group in enumerate(expression_profiles_by_labeling):
            if cur_id in patients_by_labeling[k]:
                expression_profiles_by_labeling[k].append(cur)
                patient_found = True
        if not patient_found:
            logger.info("no labeling option were found for {}".format(expression_profiles_formatted[i][0]))

    print "done split expression"

    return expression_profiles_by_labeling

def nCk(n,k):
  return long( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )

def save_sets(st,fl_name):
    fl = open(fl_name,'w+')
    for cur in st:
        for n in cur:
            fl.write(str(n)+constants.SEPARATOR)
        fl.write('\n')

def load_sets(fl_name):
    lst = []
    fl = open(fl_name,'r')
    for cur in fl.readlines():
        lst.append(cur.split(constants.SEPARATOR)[:-1])
        lst[-1][-1] = float(lst[-1][-1])
    return lst