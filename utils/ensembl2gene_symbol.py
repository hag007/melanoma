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
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
BASE_PROFILE="D:\\omics"
BASE_OUTPUT_DIR = "c:\\users\hagai\\desktop\\"


def load_gene_dictionary(gene_list_file_name, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(BASE_PROFILE,source,dataset,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines


def get_dictionary():
    lines_dict = load_gene_dictionary("ensembl2gene_symbol.txt")

    gene_symbols2ensembl = {}
    included_genes = []
    for cur in lines_dict:
        splited_line = cur.split()
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
        gene_symbols2ensembl[splited_line[1]] = splited_line[0][:limit]
    return gene_symbols2ensembl
