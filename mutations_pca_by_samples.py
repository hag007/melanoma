
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
def mutation_pca_by_samples(tested_gene_list_file_name, total_gene_list_file_name, gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, gene_filter_file_name=None, tested_gene_list_path=None, total_gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_file_path=None, var_th_index=None, is_unsupervised=True, start_k=2, end_k=2, meta_groups=None, filter_expression=None, is_ge_integ=False):

    integ_data = load_integrated_mutation_data(
                                      mutation_file_name=mutation_file_name,
                                      phenotype_file_name=phenotype_file_name,
                                      survival_file_name=survival_file_name, var_th_index=var_th_index,
                                      meta_groups=meta_groups, filter_expression=filter_expression)
    if integ_data is None:
        print "insufficient data"
        return
    mu_data, mu_data_headers_rows, mu_data_headers_columns, labels_assignment, survival_dataset, phenotype_heatmap = integ_data



    all_patients = np.unique(mu_data_headers_rows).flatten()
    all_mutated_genes = np.unique(mu_data[:,0]).flatten()
    # mis_mutated_genes = np.unique(mu_data[np.where(np.core.defchararray.find(mu_data[1:, 8], "missense")!=-1), 1]).flatten()

    all_mutated_vectors = np.zeros((len(all_patients), len(all_mutated_genes)))
    # mis_mutated_vectors = np.array([[0 for y in mis_mutated_genes] for x in range(len(all_patients))])

    print "build vectors from {} entries".format(len(mu_data))

    stopwatch = Stopwatch()
    stopwatch.start()
    a = list(all_patients)
    b = list(all_mutated_genes)
    for i, x in enumerate(mu_data):
        all_mutated_vectors[a.index(mu_data_headers_rows[i])][b.index(x[0])] += 1
    print stopwatch.stop("end mut")
    all_mutated_vectors[all_mutated_vectors>5] =5

    all_mutated_vectors = all_mutated_vectors[:,all_mutated_genes!="TTN" ]
    all_mutated_genes = all_mutated_genes[all_mutated_genes != "TTN"]

    all_mutated_vectors[all_mutated_vectors > 5] = 5
    all_mutated_genes = all_mutated_genes[(all_mutated_vectors != 0).sum(axis=0) > np.shape(all_mutated_vectors)[0] * 0.1]
    all_mutated_vectors = all_mutated_vectors[:,(all_mutated_vectors != 0).sum(axis=0) > np.shape(all_mutated_vectors)[0] * 0.1]
    print "all_mutated_vectors after filter sparse: {}".format(np.shape(all_mutated_vectors))

    if np.size(all_mutated_genes) == 0:
        return


    mutation_expression_integ = all_mutated_vectors
    mutual_patients = all_patients
    mutation_expression_integ_headers_columns = all_mutated_genes

    if is_ge_integ:
        ge_data = load_integrated_ge_data(tested_gene_list_file_name=tested_gene_list_file_name,
                                          total_gene_list_file_name=total_gene_list_file_name,
                                          gene_expression_file_name=gene_expression_file_name,
                                          phenotype_file_name=phenotype_file_name,
                                          survival_file_name=survival_file_name, var_th_index=var_th_index,
                                          meta_groups=meta_groups, filter_expression=filter_expression)
        if ge_data is None:
            print "insufficient data"
            return
        gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns, labels_assignment, survival_dataset = ge_data

        all_mutated_vectors = zscore(all_mutated_vectors, axis=0)
        gene_expression_top_var = zscore(gene_expression_top_var, axis=0)
        mutual_patients = np.array([x for x in all_patients if x in gene_expression_top_var_headers_rows])
        mutual_mutations = all_mutated_vectors[np.in1d(all_patients, mutual_patients)]
        mutual_mutations = mutual_mutations[mutual_patients.argsort()]
        mutual_patients = np.array([x for x in gene_expression_top_var_headers_rows if x in all_patients])
        mutual_expressions = gene_expression_top_var[np.in1d(gene_expression_top_var_headers_rows, mutual_patients)]
        mutual_expressions = mutual_expressions[mutual_patients.argsort()]

        mutual_patients.sort()
        mutation_expression_integ = np.c_[mutual_mutations, mutual_expressions]
        mutation_expression_integ_headers_columns = np.r_[all_mutated_genes, e2g_convertor(gene_expression_top_var_headers_columns)]
    else:
        survival_dataset = np.array(load_survival_data(survival_file_name))

        plot_pca_by_samples(mutation_expression_integ, labels_assignment, meta_groups)

