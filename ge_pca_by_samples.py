
from utils.ensembl2gene_symbol import e2g_convertor
from matplotlib import style
style.use("ggplot")
from scipy.stats import zscore
import scipy
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
import svm
from utils.km_curve import km_curve
from utils.clustering import find_clusters
from utils.pca import plot_pca_by_samples
############################ () cluster and enrichment #############################


# () main
def genes_pca_by_samples(tested_gene_list_file_name, total_gene_list_file_name, gene_expression_file_name, phenotype_file_name, survival_file_name, gene_filter_file_name=None, tested_gene_list_path=None, total_gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_file_path=None, var_th_index=None, meta_groups=None, filter_expression=None, n_components=3):
    data = load_integrated_ge_data(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, meta_groups=meta_groups, filter_expression=filter_expression)
    if data is None:
        print "insufficient data"
        return
    gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns, labels_assignment, survival_dataset = data
    plot_pca_by_samples(gene_expression_top_var, labels_assignment, meta_groups, tested_gene_list_file_name, n_components=n_components)

