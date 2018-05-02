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
from sklearn import svm
from sklearn import svm
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.model_selection import PredefinedSplit
from sklearn.metrics import accuracy_score
import scipy
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import time
from matplotlib.ticker import FormatStrFormatter
import math
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
BASE_PROFILE="D:\\omics"
BASE_OUTPUT_DIR = "c:\\users\hagai\\desktop\\"

LABEL_ID = "sample_type.samples"
PRIMARY_TUMOR = "Primary Tumor"
METASTATIC = "Metastatic"


def load_gene_list(gene_list_file_name, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(BASE_PROFILE,source,dataset,"list",gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines

def load_gene_dictionary(gene_list_file_name, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(BASE_PROFILE,source,dataset,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines

lines = load_gene_list("lysosome_genes_gene_symbols.txt")
lines_uppered = []
for i, cur in enumerate(lines):
    lines[i] = lines[i].upper()

lines_dict = load_gene_dictionary("ensembl2gene_symbol.txt")


print "genes list size: {}".format(len(lines))
print "dict genes size: {}".format(len(lines_dict))

included_genes = []
for cur in lines_dict[1:]:
    splited_line = cur.split()
    if splited_line[0].find('.') > 0:
        limit = splited_line[0].find('.')
    else:
        limit = len(splited_line[0])
    if splited_line[1] in lines:
        included_genes.append(splited_line[0][:limit])

for cur in included_genes:
    print cur

print "total:{}".format(len(included_genes))