import re
import os
import math
from numpy import *
import numpy.random
from sklearn.datasets import fetch_mldata
import sklearn.preprocessing
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
import logging
import constants
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
BASE_PROFILE="D:\\omics"
BASE_OUTPUT_DIR = "c:\\users\hagai\\desktop\\"

g2e_dict = None
e2g_dict = None
ensembl2entrez_dict = None
entrez2ensembl_dict = None

def load_gene_dictionary(gene_list_file_name, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.DICT_DIR,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines


def get_entrez2ensembl_dictionary():
    lines_dict = load_gene_dictionary(constants.ENSEMBL_TO_ENTREZ)

    entrez2ensembl = {}
    for cur in lines_dict:
        splited_line = cur.split()
        if len(splited_line) != 2: continue
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
            entrez2ensembl[splited_line[1]] = splited_line[0][:limit]
    return entrez2ensembl

def get_ensembl2entrez_dictionary():

    lines_dict = load_gene_dictionary(constants.ENSEMBL_TO_ENTREZ)
    ensembl2gene_symbols = {}
    for cur in lines_dict:
        splited_line = cur.split()
        if len(splited_line) !=2: continue
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
        ensembl2gene_symbols[splited_line[0][:limit]] = splited_line[1]
    return ensembl2gene_symbols


def ensembl2entrez_convertor(e_ids):
    global ensembl2entrez_dict
    if ensembl2entrez_dict is None:
        ensembl2entrez_dict = get_ensembl2entrez_dictionary()
    results = []
    for cur in e_ids:
        if ensembl2entrez_dict.has_key(cur.split(".")[0]):
            results.append(ensembl2entrez_dict[cur.split(".")[0]])
    return results


def entrez2ensembl_convertor(entrez_ids):
    global entrez2ensembl_dict
    if entrez2ensembl_dict is None:
        entrez2ensembl_dict = get_entrez2ensembl_dictionary()
    results = []
    for cur in entrez_ids:
        if entrez2ensembl_dict.has_key(cur.split(".")[0]):
            results.append(entrez2ensembl_dict[cur.split(".")[0]])
    return results