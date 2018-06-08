import numpy as np
import constants
import sys
from infra import *
import os

# return gene expression table filtered according an external list with proper orientation (0 degree angle according genes, 90 degree angle according patients)
def load_gene_expression_profile(gene_list_file_name=None, gene_expression_file_name=None, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None, source="GDC-TCGA",dataset="melanoma",by_gene=True, normalization_method="l_inf_norm"):
    stopwatch = Stopwatch()
    stopwatch.start()
    expression_profiles_filtered = None

    if gene_expression_path == None:
        gene_expression_path = os.path.join(constants.TCGA_DATA_DIR, gene_expression_file_name)
        stopwatch.start()
    f = open(gene_expression_path, 'r')

    if gene_list_file_name:
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

        expression_profiles_filtered = [l.strip().split() for i, l in enumerate(f) if i==0 or l[:l.strip().find('\t')].split(".")[0] in gene_list]
        # or l.strip()[0:l.strip().find('\t')] in gene_list or l.strip()[0:l.strip().find('\t')].split(".")[0] in gene_list
    else:
        f = open(gene_expression_path, 'r')
        expression_profiles_filtered = [l.strip().split() for i, l in enumerate(f)]
        # or l.strip()[0:l.strip().find('\t')] in gene_list or l.strip()[0:l.strip().find('\t')].split(".")[0] in gene_list
    print stopwatch.stop("done filter gene expression")
    f.close()
    axis = 1
    reshape_size= (len(expression_profiles_filtered) -1, 1)
    if not by_gene:
        axis= 0
        reshape_size = (1, len(expression_profiles_filtered[0])-1)
    expression_profiles_filtered = np.array(expression_profiles_filtered)
    row_header = expression_profiles_filtered[0]
    column_header = expression_profiles_filtered[1:,0]
    expression_stripped = expression_profiles_filtered[1:, 1:].astype(float)
    thismodule = sys.modules[__name__]
    norm = getattr(thismodule, normalization_method)(expression_stripped, axis)
    expression_profiles_normalized = np.nan_to_num(expression_stripped / norm.reshape(reshape_size))
    expression_profiles_normalized = np.c_[column_header, expression_profiles_normalized]
    expression_profiles_normalized = np.r_[row_header.reshape(1,len(row_header)), expression_profiles_normalized]
    return expression_profiles_normalized

def l1_norm(X, axis):
    return np.abs(X).sum(axis=axis)

def l2_norm(X, axis):
    return np.sqrt((X * X).sum(axis=axis))

def l_inf_norm(X, axis):
    return np.abs(X).max(axis=axis)

def save_expression_profile(expression_profile, output_file_name):
    f = open(os.path.join(constants.OUTPUT_DIR,output_file_name),'w+')
    for i, cur in enumerate(expression_profile):
        f.write("\t".join(cur)+"\n") # [:5]
        # if i==5000:
        #     break
    f.close()




filename = "TCGA-SKCM.htseq_fpkm{}.tsv"

is_by_gene = True

suffix = "_normalized"
if is_by_gene:
    suffix+="_by_genes"
else:
    suffix+="_by_patients"

normalized_expression_profile = load_gene_expression_profile(gene_expression_file_name=filename.format(""), by_gene=is_by_gene)
save_expression_profile(normalized_expression_profile, filename.format(suffix))