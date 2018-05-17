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
from constants import *
############################################ infra ########################################

def load_gene_list(gene_list_file_name, gene_list_path=None): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(LIST_DIR,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines

# return gene expression table filtered according an external list with proper orientation (0 degree angle according genes, 90 degree angle according patients)
def load_gene_expression_profile(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None, source="GDC-TCGA",dataset="melanoma",by_gene=False):
    stopwatch = Stopwatch()
    stopwatch.start()
    gene_list = load_gene_list(gene_list_file_name=gene_list_file_name, gene_list_path=gene_list_path)
    gene_list = [l.split(".")[0] for i, l in enumerate(gene_list)]
    print stopwatch.stop("done loading gene list")
    # random.shuffle(gene_list)
    # gene_list = gene_list[:400]
    if gene_filter_file_name:
        stopwatch.start()
        filter_gene_list = load_gene_list(gene_list_file_name=gene_filter_file_name, gene_list_path=gene_filter_path,
                                          source=source, dataset=dataset)
        gene_list = [cur for cur in gene_list if cur in filter_gene_list]
        print stopwatch.stop("done filter gene list")

    if gene_expression_path == None:
        gene_expression_path = os.path.join(TCGA_DATA_DIR, gene_expression_file_name)
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

    expression_profiles_filtered
    return expression_profiles_filtered


def load_gene_expression_profile_by_genes(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None, source="GDC-TCGA",dataset="melanoma"):
    return load_gene_expression_profile(gene_list_file_name=gene_list_file_name, gene_expression_file_name=gene_expression_file_name, gene_filter_file_name=gene_filter_file_name, gene_list_path=gene_list_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path, source=source,dataset=dataset,by_gene=True)

def load_gene_expression_profile_by_patients(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None, source="GDC-TCGA",dataset="melanoma"):
    return load_gene_expression_profile(gene_list_file_name=gene_list_file_name, gene_expression_file_name=gene_expression_file_name, gene_filter_file_name=gene_filter_file_name, gene_list_path=gene_list_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path, source=source,dataset=dataset,by_gene=False)


def load_phenotype_data(phenotype_file_name, phenotype_list_path=None, source="GDC-TCGA",dataset="melanoma"):
    if not phenotype_list_path:
        phenotype_list_path = os.path.join(TCGA_DATA_DIR,phenotype_file_name)
    f = open(phenotype_list_path, 'r')
    phenotype_profile = [l.strip().split('\t') for l in f]
    f.close()
    return phenotype_profile


def divided_patient_ids_by_tumor_type(phenotype_list_file_name, phenotype_list_path=None, source="GDC-TCGA",dataset="melanoma"):
    phenotype_data_formatted = load_phenotype_data(phenotype_list_file_name, phenotype_list_path)
    headers = phenotype_data_formatted[0]
    phenotype_profiles = phenotype_data_formatted[1:]
    primary_tumor_patients = []
    metastatic_patients = []
    label_index = [i for i, v in enumerate(headers) if v == LABEL_ID][0]
    for pp in phenotype_profiles:
        if pp[label_index] == PRIMARY_TUMOR:
            primary_tumor_patients.append(pp[0])
        elif pp[label_index] == METASTATIC:
            metastatic_patients.append(pp[0])
    return (primary_tumor_patients,metastatic_patients)

# def transform_to_map(tbl):
#     map = {}
#     for cur in tbl[1:]:
#         map[cur[0]] = cur[1:]
#     return map

# load expression profile filtered by an external genes list and divided according tumor type label
def load_expression_profile_by_gene_and_tumor_type(gene_list_file_name, gene_expression_file_name, phenotype_file_name, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_path=None, source="GDC-TCGA",dataset="melanoma"):

    expression_profiles_formatted = load_gene_expression_profile_by_patients(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=gene_filter_file_name, gene_list_path=gene_list_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path, source=source,dataset=dataset)
    primary_patients, metastatic_patients = divided_patient_ids_by_tumor_type(phenotype_file_name, phenotype_path)

    expression_profiles_primary = []
    expression_profiles_metastatic = []
    logger.info("about to split expression by primary tumor and metastatic")

    logger.info("expression_profile size: {},{}".format(*np.shape(expression_profiles_formatted)))
    for i,cur in enumerate(expression_profiles_formatted):
        if i==0: # that is, vertical headers
            expression_profiles_primary.append(cur)
            expression_profiles_metastatic.append(cur)
        elif expression_profiles_formatted[i][0] in primary_patients:
            expression_profiles_primary.append(cur)
        elif expression_profiles_formatted[i][0] in metastatic_patients:
            expression_profiles_metastatic.append(cur)
        else:
            logger.info("no tumor type for {}".format(expression_profiles_formatted[i][0]))

    print "done split expression"

    return (expression_profiles_primary, expression_profiles_metastatic)

def nCk(n,k):
  return long( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )

def save_sets(st,fl_name):
    fl = open(fl_name,'w+')
    for cur in st:
        for n in cur:
            fl.write(str(n)+SEPARATOR)
        fl.write('\n')

def load_sets(fl_name):
    lst = []
    fl = open(fl_name,'r')
    for cur in fl.readlines():
        lst.append(cur.split(SEPARATOR)[:-1])
        lst[-1][-1] = float(lst[-1][-1])
    return lst