

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
from utils.clustering import plot_heatmap
############################ () cluster and enrichment #############################


# () main

def predict_ge_by_mutation(tested_gene_list_file_name, total_gene_list_file_name, gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, gene_filter_file_name=None, tested_gene_list_path=None, total_gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_file_path=None, var_th_index=None, is_unsupervised=True, start_k=2, end_k=2, meta_groups=None, phenotype_labels_heatmap = None, filter_expression=None, integ=False, min_ratio=0.1 , included_mutation_gene_list=None, excluded_mutation_gene_list=None):


    integ_data = load_integrated_mutation_data(
                                      mutation_file_name=mutation_file_name,
                                      phenotype_file_name=phenotype_file_name,
                                      survival_file_name=survival_file_name, var_th_index=var_th_index,
                                      meta_groups=meta_groups, phenotype_labels_heatmap = phenotype_labels_heatmap,
                                      filter_expression=filter_expression)
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

    if included_mutation_gene_list is not None:
        included_mutation_gene = load_gene_list(included_mutation_gene_list)
        all_mutated_vectors = all_mutated_vectors[:, np.in1d(all_mutated_genes,included_mutation_gene)]
        all_mutated_genes = all_mutated_genes[np.in1d(all_mutated_genes,included_mutation_gene)]

    if excluded_mutation_gene_list is not None:
        excluded_mutation_gene = load_gene_list(excluded_mutation_gene_list)
        for cur in excluded_mutation_gene:
            all_mutated_vectors = all_mutated_vectors[:,all_mutated_genes!=cur ]
            all_mutated_genes = all_mutated_genes[all_mutated_genes!=cur ]

    all_mutated_vectors[all_mutated_vectors > 5] = 5
    all_mutated_genes = all_mutated_genes[(all_mutated_vectors != 0).sum(axis=0) > np.shape(all_mutated_vectors)[0] * min_ratio]
    all_mutated_vectors = all_mutated_vectors[:,(all_mutated_vectors != 0).sum(axis=0) > np.shape(all_mutated_vectors)[0] * min_ratio]
    print "all_mutated_vectors after filter sparse: {}".format(np.shape(all_mutated_vectors))

    if np.size(all_mutated_genes) == 0:
        return


    mutation_expression_integ = all_mutated_vectors
    mutual_patients = all_patients
    mutation_expression_integ_headers_columns = all_mutated_genes
    mutual_phenotype_heatmap = phenotype_heatmap
    # mutation_expression_integ = zscore(mutation_expression_integ, axis=0)
    if integ:
        ge_data = load_integrated_ge_data(tested_gene_list_file_name=tested_gene_list_file_name,
                                          total_gene_list_file_name=total_gene_list_file_name,
                                          gene_expression_file_name=gene_expression_file_name,
                                          phenotype_file_name=phenotype_file_name,
                                          survival_file_name=survival_file_name, var_th_index=var_th_index,
                                          meta_groups=meta_groups, filter_expression=filter_expression)
        if ge_data is None:
            print "insufficient data"
            return
        gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns, _ , survival_dataset = ge_data

        all_mutated_vectors = np.nan_to_num(zscore(all_mutated_vectors, axis=0))
        gene_expression_top_var = np.nan_to_num(zscore(gene_expression_top_var, axis=0))

        mutual_patients = np.array([x for x in all_patients if x in gene_expression_top_var_headers_rows])
        mutual_mutations = all_mutated_vectors[np.in1d(all_patients, mutual_patients)]
        mutual_mutations = mutual_mutations[mutual_patients.argsort()]

        if phenotype_labels_heatmap is not None:
            mutual_phenotype_heatmap = phenotype_heatmap[np.in1d(all_patients, mutual_patients)]
            mutual_phenotype_heatmap = mutual_phenotype_heatmap[mutual_patients.argsort()]

        mutual_patients = np.array([x for x in gene_expression_top_var_headers_rows if x in all_patients])
        mutual_expressions = gene_expression_top_var[np.in1d(gene_expression_top_var_headers_rows, mutual_patients)]
        mutual_expressions = mutual_expressions[mutual_patients.argsort()]

        mutual_patients.sort()
        mutation_expression_integ = np.c_[mutual_mutations, mutual_expressions]
        mutation_expression_integ_headers_columns = np.r_[all_mutated_genes, gene_expression_top_var_headers_columns]
    else:
        survival_dataset = np.array(load_survival_data(survival_file_name))

    if is_unsupervised:
        print "find clusters"
        clfs_results = find_clusters(end_k, mutation_expression_integ, mutual_patients,
                                    start_k, mutation_expression_integ_headers_columns,
                                   tested_gene_list_file_name, labels_assignment=labels_assignment, phenotype_heatmap=mutual_phenotype_heatmap)
        for cur_k in range(start_k, end_k+1):
            km_curve(clfs_results[cur_k], survival_dataset[1:], mutual_patients,
                     tested_gene_list_file_name.split(".")[0], i)
    else:
        for i, cur_groups in enumerate(meta_groups):
            labeled_patients = divided_patient_ids_by_label(phenotype_file_name, groups=cur_groups)
            plot_heatmap(mutation_expression_integ, mutation_expression_integ_headers_columns,
                         [labels_assignment[i]] + labels_assignment[:i] + labels_assignment[i + 1:],
                         mutual_patients,
                         tested_gene_list_file_name, label_index=i)
            km_curve(labeled_patients, survival_dataset[1:], mutual_patients, tested_gene_list_file_name.split(".")[0],label_index=i)



    # print "assign labels"
    # labels = np.array([0 for x in all_patients])
    # for i, cur_group in enumerate(clfs_results[2]):
    #     labels[np.in1d(all_patients, cur_group)] = i
    # r_prediction_mutation, r_prediction_labels = svm.randonize_patients(all_mutated_vectors, labels)
    #
    #
    # data_train, data_test, labels_train, labels_test = svm.divide_train_and_test_groups(r_prediction_mutation, r_prediction_labels)
    # tuning_parameters = {'C': [10], 'kernel': ['rbf']}
    # clf = svm.svm_rbf_default(tuning_parameters)
    #
    # for cur in range(10):
    #     r_prediction_mutation, r_prediction_labels = svm.randonize_patients(all_mutated_vectors, labels)
    #     data_train, data_test, labels_train, labels_test = svm.divide_train_and_test_groups(r_prediction_mutation, r_prediction_labels)
    #     print "cur iteration: {}".format(cur)
    #     test_pr, test_roc = svm.apply_svm(clf, data_train, labels_train, data_test, labels_test, "DISTANCE")
    #
    # print labels


