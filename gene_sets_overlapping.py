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
def find_sets_overlaps(tested_gene_list_file_name, total_gene_list_file_name, gene_expression_file_name, phenotype_file_name, gene_filter_file_name=None, tested_gene_list_path=None, total_gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_file_path=None, source="GDC-TCGA",dataset="melanoma", hgt_preprocessing_file_name = None, pval_preprocessing_file_name= None, N = None, B = None):
    print "about ot analyse: {}".format(tested_gene_list_file_name)
    # fetch gene expression by gene_id, divided by tumor type
    gene_sets = []
    for cur in tested_gene_list_file_name:
        gene_sets.append(load_gene_list(cur))

    output = ""
    for i, cur_1 in enumerate(gene_sets):
        output += "\tgroup {} (n={}):".format(tested_gene_list_file_name[i].split('.')[0], len(gene_sets[i]))
    output += "\n"
    for i, cur_1 in enumerate(gene_sets):
        output += "group {} (n={}):".format(tested_gene_list_file_name[i].split('.')[0], len(gene_sets[i]))
        for j, cur_2 in enumerate(gene_sets):
            if j <= i:
                output += "\t"+str(len([x for x in cur_1 if x in cur_2]))
        output += "\n"

    f = file(os.path.join(constants.OUTPUT_DIR, "OVERLAPPING_GENE_SETS_{}_{}.txt".format(gene_expression_file_name.split(".")[0], time.time())), 'w+')
    f.write(output)
    f.close()
