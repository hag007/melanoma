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

class GeneAnalyzer:

    ############################################ infra ########################################
    def __init__(self, tested_gene_list_file_name, gene_expression_file_name, phenotype_file_name, filter_gene_list_file_name=None,tested_gene_list_path=None, gene_expression_path=None, phenotype_path=None, filter_gene_list_path=None, source="GDC-TCGA",dataset="melanoma"):
        self.tested_gene_list_file_name = tested_gene_list_file_name
        self.gene_expression_file_name=gene_expression_file_name
        self.phenotype_file_name=phenotype_file_name
        self.filter_gene_list_file_name=filter_gene_list_file_name
        self.tested_gene_list_path=tested_gene_list_path
        self.gene_expression_path=gene_expression_path
        self.phenotype_path=phenotype_path
        self.filter_gene_list_path=filter_gene_list_path
        self.source="GDC-TCGA"
        self.dataset="melanoma"

    def load_gene_list(self, gene_list_file_name=None, gene_list_path=None, source=None,dataset=None):
        gene_list_file_name = self.gene_list_file_name, gene_list_path = self.gene_list_path, source = self.source, dataset = self.dataset

        if gene_list_path == None:
            gene_list_path = os.path.join(BASE_PROFILE,source,dataset,"list",gene_list_file_name)
        f = open(gene_list_path,'r')
        lines = [l.strip() for l in f]
        f.close()
        return lines

    # return gene expression table filtered according an external list with proper orientation (0 degree angle according genes, 90 degree angle according patients)
    def load_gene_expression_profile(self, gene_list_file_name=self.gene_list_file_name, gene_expression_file_name=self.gene_expression_file_name,
                                     gene_list_path=self.gene_list_path, gene_expression_path=self.gene_expression_path,
                                     source=self.source,dataset=self.dataset,by_gene=False, gene_filter_file_name=self.gene_filter_file_name):

        gene_list = load_gene_list()
        if gene_filter_file_name:
            filter_gene_list = load_gene_list(gene_list_file_name=gene_filter_file_name, gene_list_path=gene_filter_path)
            gene_list = [cur for cur in gene_list if cur in filter_gene_list or cur[:cur.find('.')] in filter_gene_list]

        if gene_expression_path == None:
            gene_expression_path = os.path.join(BASE_PROFILE, source, dataset, gene_expression_file_name)
            f = open(gene_expression_path,'r')
        expression_profile_raw = [l.strip().split() for i, l in enumerate(f) if i==0 or any([l.strip()[0:l.strip().find('\t')] in gene_list or l.strip()[0:l.strip().find('\t')].split(".")[0] in gene_list])]
        f.close()
        if not by_gene:
            expression_profiles_filtered = np.flip(np.rot90(expression_profile_raw, k=1, axes=(1,0)),1)

        return expression_profiles_filtered


    def load_gene_expression_profile_by_genes(self, gene_list_file_name=self.gene_list_file_name, gene_expression_file_name=self.gene_expression_file_name, gene_filter_file_name=self.gene_filter_file_name,
                                              gene_list_path=self.gene_list_path, gene_expression_path=self.gene_expression_path, gene_filter_path=self.gene_filter_path,
                                     source=self.source,dataset=self.dataset):
        return load_gene_expression_profile(gene_list_file_name=gene_list_file_name, gene_expression_file_name=gene_expression_file_name, gene_filter_file_name=gene_filter_file_name,
                                              gene_list_path=gene_list_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path,
                                     source=source,dataset=dataset, by_genes=True)


    def load_gene_expression_profile_by_patients(self, gene_list_file_name=self.gene_list_file_name, gene_expression_file_name=self.gene_expression_file_name, gene_filter_file_name=self.gene_filter_file_name,
                                              gene_list_path=self.gene_list_path, gene_expression_path=self.gene_expression_path, gene_filter_path=self.gene_filter_path,
                                     source=self.source,dataset=self.dataset):
        return load_gene_expression_profile(gene_list_file_name=gene_list_file_name, gene_expression_file_name=gene_expression_file_name, gene_filter_file_name=gene_filter_file_name,
                                              gene_list_path=gene_list_path, gene_expression_path=gene_expression_path, gene_filter_path=gene_filter_path,
                                     source=source,dataset=dataset, by_genes=False)


    def load_phenotype_data(self, phenotype_file_name=self.phenotype_file_name, phenotype_list_path=self.phenotype_list_path, source=self.source,dataset=self.dataset):
        phenotype_list_path = os.path.join(BASE_PROFILE,source,dataset,phenotype_file_name)
        f = open(phenotype_list_path, 'r')
        phenotype_profile = [l.strip().split('\t') for l in f]
        f.close()
        return phenotype_profile


    def divided_patient_ids_by_tumor_type(self, phenotype_list_file_name=self.phenotype_list_file_name, phenotype_list_path=self.phenotype_list_path, source=self.source,dataset=self.dataset):
        phenotype_data_formatted = load_phenotype_data(phenotype_list_file_name)
        headers = phenotype_data_formatted[0]
        phenotype_profiles = phenotype_data_formatted[1:]
        primary_tumor_patients = []
        metastatic_patients = []
        label_index = [i for i, v in enumerate(headers) if v == LABEL_ID][0]
        for pp in phenotype_profiles:
            if pp[label_index] == PRIMARY_TUMOR:
                primary_tumor_patients.append(pp[0])
            elif pp[label_index] == METASTATIC:
                metastatic_patients.append(pp[0])
        return (primary_tumor_patients,metastatic_patients)

    # def transform_to_map(tbl):
    #     map = {}
    #     for cur in tbl[1:]:
    #         map[cur[0]] = cur[1:]
    #     return map

    # load expression profile filtered by an external genes list and divided according tumor type label
    def load_expression_profile_by_gene_and_tumor_type(self, gene_list_file_name=self.gene_list_file_name, gene_expression_file_name=self.gene_expression_file_name, phenotype_file_name=self.phenotype_file_name, gene_list_path=self.gene_list_path, gene_expression_path=self.gene_expression_path, phenotype_path=self.phenotype_path, source=self.source, dataset=self.dataset, gene_filter_file_name=self.gene_filter_file_name):

        expression_profiles_formatted = load_gene_expression_profile_by_patients(gene_list_file_name, gene_expression_file_name, gene_filter_file_name=None, gene_list_path=None, gene_expression_path=None, gene_filter_path=None, source="GDC-TCGA",dataset="melanoma")
        primary_patients, metastatic_patients = divided_patient_ids_by_tumor_type(phenotype_file_name)

        expression_profiles_primary = []
        expression_profiles_metastatic = []
        logger.info("about to split expression by primary tumor and metastatic")

        logger.info("expression_profile size: {},{}".format(*np.shape(expression_profiles_formatted)))
        for i,cur in enumerate(expression_profiles_formatted):
            if i==0: # that is, vertical headers
                expression_profiles_primary.append(cur)
                expression_profiles_metastatic.append(cur)
            elif expression_profiles_formatted[i][0] in primary_patients:
                expression_profiles_primary.append(cur)
            elif expression_profiles_formatted[i][0] in metastatic_patients:
                expression_profiles_metastatic.append(cur)
            else:
                logger.info("no tumor type for {}".format(expression_profiles_formatted[i][0]))

        print "done split expression"

        return (expression_profiles_primary, expression_profiles_metastatic)


    ################################################# (1) eucalidian distance between patients according tested genes expression profile #############################################

    # (1) main
    def find_expression_similarity_profile(self):
        primary_labeled, metastatic_labeled = load_expression_profile_by_gene_and_tumor_type()

        similarity_primary = []
        similarity_metastatic = []
        inter_similarity = []
        counter = 0
        for cur_1 in primary_labeled[1:]:
            counter += 1
            for cur_2 in primary_labeled[1:]:
                similarity_primary.append(np.linalg.norm(np.array(cur_1[1:], dtype=float)-np.array(cur_2[1:], dtype=float)))

        counter=0
        for cur_1 in metastatic_labeled[1:]:
            counter+=1
            for cur_2 in metastatic_labeled[1:]:
                similarity_metastatic.append(np.linalg.norm(np.array(cur_1[1:], dtype=float)-np.array(cur_2[1:], dtype=float)))

        counter = 0
        for cur_1 in primary_labeled[1:]:
            counter += 1
            for cur_2 in metastatic_labeled[1:]:
                inter_similarity.append(np.linalg.norm(np.array(cur_1[1:], dtype=float)-np.array(cur_2[1:], dtype=float)))

        print "similarity_primary: {}, {}".format(np.average(similarity_primary), np.var(similarity_primary))
        print "similarity_metastatic: {}, {}".format(np.average(similarity_metastatic), np.var(similarity_metastatic))
        print "inter_similarity: {}, {}".format(np.average(inter_similarity), np.var(inter_similarity))
        return np.average(similarity_primary), np.var(similarity_primary), np.average(similarity_metastatic), np.var(similarity_metastatic), np.average(inter_similarity), np.var(inter_similarity)

    ############################ (2) significance expression and proportion differntiations #############################


    def plot_pvalues(self, y_axis, x_axis, th,output_file_name, is_bigger_better=False):
        n, bins, patches = mpl.pyplot.hist(y_axis, x_axis)
        for c, p in zip(bins, patches):
            if is_bigger_better:
                th_condition = c < th
            else:
                th_condition = c > th

            if th_condition:
                color = 'blue'
            else:
                color = 'red'
            plt.setp(p, 'facecolor', color)
        plt.savefig(os.path.join(BASE_OUTPUT_DIR, output_file_name))
        plt.cla()


    def plot_pvalues_log_scaled(self, y_axis, x_axis, th,output_file_name):
        plot_pvalues([-math.log(cur, 10) for cur in y_axis], x_axis, -math.log(th, 10), output_file_name, is_bigger_better=True)


    def summarize_genes_proportion(tested_gene_list, total_gene_list, gene_pval_pair, true_counter):
        significant_tested_gene_list = [((i + 1), cur[0], cur[1]) for i, cur in enumerate(gene_pval_pair) if
                                        cur[0] in tested_gene_list and i < true_counter]
        included_tested_genes = [cur for cur in tested_gene_list if
                                 cur in total_gene_list or cur[:cur.find('.')] in total_gene_list]
        included_tested_genes_size = len(included_tested_genes)
        significant_tested_gene_list_size = len(significant_tested_gene_list)
        print "total tested genes in true hypothsis: {} out of possible {}".format(significant_tested_gene_list_size, included_tested_genes_size)
        tested_gene_list_size = len(tested_gene_list)
        total_gene_list_size = len(total_gene_list)
        results_table = []
        expected_actual_difference_list = []
        expected_actual_ratio_difference_list = []
        rank_in_total_list_list = []
        z_test_proportion_test_list = []
        for i, cur in enumerate(significant_tested_gene_list):
            rank_in_tesed_list = (i + 1)
            rank_in_total_list = cur[0]
            ensembel_id = cur[1]
            p_val = cur[2]
            expected_quantity = rank_in_total_list * (included_tested_genes_size / (total_gene_list_size * 1.0))
            expected_proportion = included_tested_genes_size / (total_gene_list_size * 1.0)
            actual_quantity = (i + 1)
            actual_proportion = (i + 1) / (cur[0] * 1.0)
            expected_actual_difference = actual_quantity - expected_quantity
            expected_actual_ratio_difference = expected_actual_difference / (expected_quantity * 1.0)
            z_test_proportion_test = (actual_proportion - expected_proportion) / math.sqrt(
                (expected_proportion * (1 - expected_proportion)) / rank_in_total_list)
            results_table.append([rank_in_tesed_list, rank_in_total_list, ensembel_id, p_val,
                                  expected_quantity, expected_proportion,
                                  actual_quantity, actual_proportion,
                                  expected_actual_difference, expected_actual_ratio_difference,
                                  z_test_proportion_test])

            expected_actual_difference_list.append(expected_actual_difference)
            expected_actual_ratio_difference_list.append(expected_actual_ratio_difference)
            z_test_proportion_test_list.append(z_test_proportion_test)
            rank_in_total_list_list.append(rank_in_total_list)
        return expected_actual_difference_list, expected_actual_ratio_difference_list, rank_in_total_list_list, z_test_proportion_test_list


    def plot_genes_proportion(self, expected_actual_difference_list, expected_actual_ratio_difference_list, z_test_proportion_test_list, rank_in_total_list_list, total_significant_hypotheses_size, expected_tested_genes_ratio, tested_gene_list_file_name):
        z_score_threshold_two_way = 1.96
        tested_genes_size = len(rank_in_total_list_list)
        y_counter = [min(i, tested_genes_size) for i in range(1,tested_genes_size+1)]
        plt.plot(rank_in_total_list_list[1:], y_counter[1:], label="number of significant values (n)")
        plt.plot(rank_in_total_list_list[1:], expected_actual_difference_list[1:], label="actual-expected significant hypo. difference (n)")
        plt.plot([total_significant_hypotheses_size, total_significant_hypotheses_size], [-20, tested_genes_size + 5], label="True hypotheses threshold")
        plt.plot([0, total_significant_hypotheses_size], [0, tested_genes_size],
                 color="gray")
        plt.plot([0, total_significant_hypotheses_size], [0, expected_tested_genes_ratio],
                 color="black")
        plt.legend()
        plt.savefig(os.path.join(BASE_OUTPUT_DIR, "{}_sum_n".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')])))
        plt.cla()
        plt.plot(rank_in_total_list_list[1:], expected_actual_ratio_difference_list[1:], label="actual/expected proportion ratio")
        plt.plot(rank_in_total_list_list[1:], z_test_proportion_test_list[1:], label="z_score")
        plt.plot(rank_in_total_list_list[1:], [z_score_threshold_two_way for i in range(1,tested_genes_size)], label="z_score threshold (two-way)")
        plt.plot([total_significant_hypotheses_size, total_significant_hypotheses_size], [-0.5, 3.5], label="True hypotheses threshold")
        plt.legend()
        plt.savefig(os.path.join(BASE_OUTPUT_DIR,"{}_sum_p".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')])))

    # (2) main
    def find_expression_significance(self, tested_gene_list_file_name, total_gene_list_file_name):

        # fetch gene expression by gene_id, divided by tumor type
        primary_expression, metastatic_expression = load_expression_profile_by_gene_and_tumor_type(total_gene_list_file_name)
        primary_expression = np.rot90(np.flip(primary_expression, 1), k=-1, axes=(1,0))
        metastatic_expression = np.rot90(np.flip(metastatic_expression, 1), k=-1, axes=(1, 0))

        # test pval for significance differentiation between label values (primar vs metastatic)
        pvals = []
        gene_symbols = []
        for  i in range(1,len(primary_expression)):
            cur_pval = scipy.stats.ttest_ind([float(c) for c in primary_expression[i][1:]], [float(c) for c in metastatic_expression[i][1: ]])[1]
            if not math.isnan(cur_pval):
                pvals.append(cur_pval)
                gene_symbols.append(primary_expression[i][0])

        # sort gene_id-pval pairs by pval
        gene_pval_pair = zip(gene_symbols,pvals)
        gene_pval_pair.sort(key=lambda x: x[1], reverse=False)

        # plot pvals
        pvals = [cur[1] for cur in gene_pval_pair]
        logger.info("pvals (uncorrected) below 0.05: {}".format(np.sum([True if pv < 0.05 else False for pv in pvals])))
        plot_pvalues(pvals, [i * 0.01 for i in range(101)], 0.05, "{}_01".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')]))
        plot_pvalues_log_scaled(pvals, range(100), 0.05, "{}_02".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')]))

        # correct pval (i.e generate qvals out of pvals)
        fdr_results = fdrcorrection0(pvals, alpha=0.05, method='indep', is_sorted=True)
        true_counter = len([cur for cur in fdr_results[0] if cur == True])
        print "true hypothesis: {}".format(true_counter)
        print "total hypothesis: {}".format(np.size(fdr_results[0]))

        # plot qvals
        qvals = fdr_results[1]
        qvals.sort()
        qval_threshold = qvals[true_counter-1]
        plot_pvalues(qvals, [i * 0.01 for i in range(101)], qval_threshold, "{}_11".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')]))
        plot_pvalues_log_scaled(qvals, range(100), qval_threshold, "{}_22".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')]))

        # n, bins, patches = mpl.pyplot.hist([cur[0] for cur in significant_oxidative_gene_list], range(0,len(total_gene_list)), cumulative=True)
        # mpl.pyplot.savefig("c:\\users\hagai\\desktop\\figure_3")
        # mpl.pyplot.cla()

        # summarize genes proportion
        tested_gene_list = load_gene_list(gene_list_file_name=tested_gene_list_file_name)
        total_gene_list = load_gene_list(gene_list_file_name=total_gene_list_file_name)
        expected_actual_difference_list, expected_actual_ratio_difference_list, rank_in_total_list_list, \
        z_test_proportion_test_list = summarize_genes_proportion(tested_gene_list, total_gene_list, gene_pval_pair, true_counter)


        # calculate_gene_concentration
        expected_tested_genes_ratio = len(tested_gene_list)*true_counter/(len(total_gene_list)*1.0)
        plot_genes_proportion(expected_actual_difference_list, expected_actual_ratio_difference_list, z_test_proportion_test_list, rank_in_total_list_list, true_counter, expected_tested_genes_ratio, tested_gene_list_file_name)

    ############################################(3) SVM to see whether we can predict tumor type by tested genes' profile expression########################################

    def load_svm_data(self, tested_gene_file_name):
        data = []
        labels = []
        primary_labeled, metastatic_labeled = load_expression_profile_by_gene_and_tumor_type(tested_gene_file_name)

        for cur in primary_labeled[1:]:
            data.append(cur[1:])
            labels.append(0)

        for cur in metastatic_labeled[1:]:
            data.append(cur[1:])
            labels.append(1)

        return data, labels

    def randonize_patients(self, data, labels,data_alt,labels_alt):
        c = list(zip(data, labels, data_alt, labels_alt))
        random.shuffle(c)
        return zip(*c)

    def divide_train_and_test_groups(self, data, labels):
        data_train = data[:(3 * len(data)) / 4]
        data_test = data[(3 * len(data)) / 4:]
        labels_train = labels[:(3 * len(data)) / 4]
        labels_test = labels[(3 * len(data)) / 4:]
        return data_train, data_test, labels_train, labels_test

    def apply_svm(self, tuned_parameters, data_train, labels_train, data_test, labels_test):
        gs_train = GridSearchCV(svm.SVC(), param_grid=tuned_parameters, return_train_score=True)
        gs_train.fit(data_train, labels_train)
        predicted_results = gs_train.predict(data_test)
        ta = (len(labels_test) * 1.0 - sum([abs(p - r) for p, r in zip(predicted_results, labels_test)])) / (len(labels_test) * 1.0)
        return gs_train.best_score_, ta

    def print_svm_results(self, train_scores, train_alt_scores, test_scores, test_alt_scores, rounds):
        print "train_accuracy_avg: {}".format(sum(train_scores)/ rounds)
        print "train_alt_accuracy_avg: {}".format(sum(train_alt_scores)/ rounds)
        print "train_accuracy_diff_avg: {}".format((sum(train_scores)-sum(train_alt_scores))/rounds)
        print "test_accuracy_avg: {}".format(sum(test_scores)/ rounds)
        print "test_alt_accuracy_avg: {}".format(sum(test_alt_scores)/ rounds)
        print "test_accuracy_diff_avg: {}".format((sum(test_scores)-sum(test_alt_scores))/rounds)
        print "train p val: {}".format(scipy.stats.ttest_ind(train_scores, train_alt_scores)[1])
        print "test_p_val: {}".format(scipy.stats.ttest_ind(test_scores, test_alt_scores)[1])

        # print "{}".format(train_score_sum / rounds)
        # print "{}".format(train_alt_score_sum / rounds)
        # print "{}".format(scipy.stats.ttest_ind(train_scores, train_alt_scores)[1])
        # print "{}".format(test_score_sum / rounds)
        # print "{}".format(test_alt_score_sum / rounds)
        # print "{}".format(scipy.stats.ttest_ind(test_scores, train_alt_scores)[1])


    # (3) main
    def predict_tumor_type_by_tested_gene_expression(self, tested_gene_file_name, tested_gene_file_name_alt, rounds=100):

        data, labels = load_svm_data(tested_gene_file_name)
        data_alt, labels_alt = load_svm_data(tested_gene_file_name_alt)

        train_scores = []
        train_alt_scores = []
        test_scores = []
        test_alt_scores = []
        for i in range(rounds):
            data, labels, data_alt, labels_alt = randonize_patients(data, labels, data_alt,labels_alt)
            data_train, data_test, labels_train, labels_test = divide_train_and_test_groups(data, labels)
            data_alt_train, data_alt_test, labels_alt_train, labels_alt_test = divide_train_and_test_groups(data_alt,labels_alt)

            tuned_parameters = {'C': [10], 'kernel': ['rbf']}
            train_accuracy, test_accuracy = apply_svm(tuned_parameters, data_train, labels_train, data_test, labels_test)
            alt_train_accuracy, alt_test_accuracy = apply_svm(tuned_parameters, data_alt_train, labels_alt_train, data_alt_test, labels_alt_test)
            # print "train_accuracy: {}, alt: {}, diff {} ".format(train_accuracy, alt_train_accuracy, alt_train_accuracy - train_accuracy)
            # print "test_accuracy: {}, alt: {}, diff {} ".format(test_accuracy, alt_test_accuracy, test_accuracy - alt_test_accuracy)

            train_scores.append(train_accuracy)
            train_alt_scores.append(alt_train_accuracy)
            test_scores.append(test_accuracy)
            test_alt_scores.append(alt_test_accuracy)
        print "#######################################"
        print "RESULTS FOR {} VS {}:".format(
            tested_gene_file_name[:tested_gene_file_name.find('.')], tested_gene_file_name_alt[:tested_gene_file_name_alt.find('.')])
        print_svm_results(train_scores, train_alt_scores, test_scores, test_alt_scores, float(rounds))

obj = GeneAnalyzer()
# (1) similarity caller:
# find_expression_similarity_profile(gene_list_file_name="oxidative.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", phenotpe_file_name="TCGA-SKCM.GDC_phenotype.tsv")
# (2) call:
obj.find_expression_significance("oxidative.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
# find_expression_significance("apoptosis_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
# find_expression_significance("lysosome_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
# find_expression_significance("lipid_synthase_genes.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
# find_expression_significance("mito.txt", "protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
# (3) call
predict_tumor_type_by_tested_gene_expression("mito.txt", "oxidative.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rounds=100)
predict_tumor_type_by_tested_gene_expression("mito.txt", "top420s.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rounds=100)
predict_tumor_type_by_tested_gene_expression("mito.txt", "apoptosis_genes.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rounds=100)
predict_tumor_type_by_tested_gene_expression("mito.txt", "lysosome_genes.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rounds=100)
predict_tumor_type_by_tested_gene_expression("mito.txt", "lipid_synthase_genes.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", rounds=100)


# primaries_avg = []
# primaries_var = []
# metstatics_avg = []
# metstatics_var = []
# inters_avg = []
# inters_var = []
#
# for i in range(100):
#     print "interation {}".format(i)
#     primary_similarity_avg, primary_similarity_var, metastatic_similarity_avg, metastatic_similarity_var, inter_similarity_avg, inter_similarity_var = find_expression_similarity_profile()
#     print "primary:{}, {}; metastatic:{}, {}; inter: {}, {}".format(primary_similarity_avg, primary_similarity_var, metastatic_similarity_avg, metastatic_similarity_var, inter_similarity_avg, inter_similarity_var)
#
#
#     primaries_avg.append(primary_similarity_avg)
#     primaries_var.append(primary_similarity_var)
#     metstatics_avg.append(metastatic_similarity_avg)
#     metstatics_var.append(metastatic_similarity_var)
#     inters_avg.append(inter_similarity_avg)
#     inters_var.append(inter_similarity_var)
#
# print primaries_avg
# print primaries_var
# print metstatics_avg
# print metstatics_var
# print inters_avg
# print inters_var


