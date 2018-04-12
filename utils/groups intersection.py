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


lines_1 = load_gene_list("lysosome_genes.txt")
lines_uppered = []
for i, cur in enumerate(lines_1):
    lines_1[i] = lines_1[i].upper()

lines_2 = load_gene_list("protein_coding.txt")
lines_uppered = []
for i, cur in enumerate(lines_2):
    lines_2[i] = lines_2[i].upper()

print "genes_1 list size: {}".format(len(lines_1))
print "genes_2 list size: {}".format(len(lines_2))

included_genes = []
for cur in lines_1:
    if cur in lines_2 or cur[0:cur.find('.')] in lines_2:
        included_genes.append(cur)

for cur in included_genes:
    print cur

print "gene count group 1: {} gene count groups two: {}, mutual:{}".format(len(lines_1), len(lines_2),len(included_genes))