import scipy.special
import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
import scipy
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from scipy.stats import pearsonr
############################ (2) correlations #############################


# () main
def find_sets_correlations(tested_gene_list_file_name, total_gene_list_file_name, gene_expression_file_name, phenotype_file_name, gene_filter_file_name=None, tested_gene_list_path=None, total_gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_file_path=None, source="GDC-TCGA",dataset="melanoma", hgt_preprocessing_file_name = None, pval_preprocessing_file_name= None, N = None, B = None):
    print "about ot analyse: {}".format(tested_gene_list_file_name)
    # fetch gene expression by gene_id, divided by tumor type
    gene_sets = []
    expression_sets = []
    averaged_expression_sets = []
    total_gene_expression = load_gene_expression_profile_by_genes(total_gene_list_file_name, gene_expression_file_name, gene_filter_file_name, total_gene_list_path, gene_expression_path, phenotype_path, gene_filter_file_path)
    for cur in tested_gene_list_file_name:
        gene_sets.append(load_gene_list(cur))
        expression_sets.append([x[1:] for x in total_gene_expression if x[0].split('.')[0] in gene_sets[-1] or x[0] in gene_sets[-1]])
        averaged_expression_sets.append((np.mean(np.array(expression_sets[-1]).astype('float32'),axis=0), np.mean(np.var(np.array(expression_sets[-1]).astype('float32'), axis=0))))

    print "\tavg var",
    for i, cur_1 in enumerate(averaged_expression_sets):
        print "\tgroup {} (n={}):".format(tested_gene_list_file_name[i].split('.')[0], len(gene_sets[i])),
    print ""
    for i, cur_1 in enumerate(averaged_expression_sets):
        print "group {} (n={}):".format(tested_gene_list_file_name[i].split('.')[0], len(gene_sets[i])),
        print "\t"+str(cur_1[1]),
        for j, cur_2 in enumerate(averaged_expression_sets):
            if j <= i:
                print "\t"+str(pearsonr(cur_1[0], cur_2[0])[0]),
        print ""

