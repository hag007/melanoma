import Bio.UniProt.GOA as GOA
# from orangecontrib.bio.go import Ontology
import wget
from utils.ensembl2gene_symbol import e2g_convertor
import time
import requests
import scipy.special
import matplotlib.pyplot as plt
from matplotlib import style
import matplotlib.ticker as ticker
style.use("ggplot")
import scipy
from lifelines.statistics import logrank_test
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
import openpyxl
from openpyxl import Workbook
from openpyxl.styles import Color, PatternFill, Font, Border, Side, Alignment
import sys
import os
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_ncbi_gene2go
# from goatools.associations import read_gaf
from lifelines import KaplanMeierFitter
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import rankdata
from genes_statistic import plot_genes_statistic

############################ () cluster and enrichment #############################


# () main
def find_clusters_and_survival(tested_gene_list_file_name, total_gene_list_file_name, gene_expression_file_name, phenotype_file_name, survival_file_name, gene_filter_file_name=None, tested_gene_list_path=None, total_gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_file_path=None, var_th_index=None, is_unsupervised=True, start_k=2, end_k=2, meta_groups=None, filter_expression=None):

    data = load_integrated_ge_data(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name,                                                                                                                                                    phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, meta_groups=meta_groups, filter_expression=filter_expression)
    if data is None:
        print "insufficient data"
        return
    gene_expression_top_var, gene_expression_top_var_headers_rows, gene_expression_top_var_headers_columns, labels_assignment, survival_dataset = data

    plot_genes_statistic(gene_expression_top_var, gene_expression_top_var_headers_columns, tested_gene_list_file_name)

    if is_unsupervised:
        clfs_results = find_clusters(end_k, gene_expression_top_var, gene_expression_top_var_headers_rows,
                                    start_k, gene_expression_top_var_headers_columns,
                                   tested_gene_list_file_name, labels_assignment)

        for i in range(start_k,end_k+1):
            km_curve(clfs_results[i], survival_dataset[1:], gene_expression_top_var_headers_rows, tested_gene_list_file_name.split(".")[0],i)
    else:
        for i, cur_groups in enumerate(meta_groups):
            labeled_patients = divided_patient_ids_by_label(phenotype_file_name, groups=cur_groups)
            km_curve(labeled_patients, survival_dataset[1:], gene_expression_top_var_headers_rows, tested_gene_list_file_name.split(".")[0],label_index=i)
            plot_heatmap(gene_expression_top_var, gene_expression_top_var_headers_rows,
                         [labels_assignment[i]] + labels_assignment[:i] + labels_assignment[i+1:], gene_expression_top_var_headers_columns,
                             tested_gene_list_file_name, label_index=i)


def find_clusters(end_k, gene_expression_top_var, gene_expression_top_var_headers_rows, start_k,
                gene_expression_top_var_headers_columns, tested_gene_list_file_name, labels_assignment=None):
    clfs_results = {}
    for n_clusters in range(start_k, end_k + 1):
        clfs_results[n_clusters] = []

        km_clf = KMeans(n_clusters).fit(gene_expression_top_var)
        ranks = []
        for i in range(n_clusters):
            ranks.append(np.average(np.delete(gene_expression_top_var, np.where(km_clf.labels_ != i)[0], axis=0)))
        ranks = rankdata(ranks)
        cluster_labels = np.array(km_clf.labels_)
        for i in range(n_clusters):
            cluster_labels[np.where(km_clf.labels_ == ranks[i] - 1)] = i
        if labels_assignment is not None:
            labels_assignment = [cluster_labels + 1] + labels_assignment
        else:
            labels_assignment = [cluster_labels + 1]

        for i in range(n_clusters):
            cluster_indices = np.where(cluster_labels != i)[0]
            gene_expression_cluster = np.delete(gene_expression_top_var_headers_rows, cluster_indices, axis=0)
            gene_headers_column_cluster = np.delete(gene_expression_top_var_headers_rows, cluster_indices, axis=0)
            clfs_results[n_clusters].append(gene_headers_column_cluster)
            desc = "k={} clustering cluster {} has {} genes".format(n_clusters, i, len(gene_expression_cluster))

        # gene_expression_top_var = np.c_[gene_expression_top_var,km_clf.labels_]

        plot_heatmap(gene_expression_top_var, gene_expression_top_var_headers_rows,
                     labels_assignment, gene_expression_top_var_headers_rows,
                     tested_gene_list_file_name, n_clusters)
    return clfs_results


def labels_assignments(meta_groups, phenotype_file_name, patients_list):
    labels_assignment = []
    for i, cur_groups in enumerate(meta_groups):
        labeled_patients = divided_patient_ids_by_label(phenotype_file_name, groups=cur_groups)
        cur_labeled = np.array(patients_list)
        for j, cur_patients_group in enumerate(labeled_patients):
            cur_labeled[np.in1d(cur_labeled, cur_patients_group)] = j + 1
        cur_labeled[~np.core.defchararray.isdigit(cur_labeled)] = 0
        labels_assignment.append(cur_labeled.astype(np.int32))
    return labels_assignment


def plot_heatmap(gene_expression_top_var, gene_expression_top_var_headers_rows, labels_assignment,
                 tested_gene_expression_headers_columns, tested_gene_list_file_name, n_clusters = None, label_index=None):
    cluster_labels = labels_assignment[0]
    fig = plt.figure()
    abs_left = 0.05
    abs_bottom = 0.05
    abs_height = 0.9
    main_width = 0.8
    label_width = 0.025
    axes = []
    ax1 = fig.add_axes([abs_left, abs_bottom, main_width, abs_height])
    axes.append(ax1)
    for i, cur in enumerate(labels_assignment):
        ax2 = fig.add_axes(
            [abs_left + main_width + label_width + (i * 2) * label_width, abs_bottom, label_width, abs_height])
        axes.append(ax2)
        ax2.grid(False)
        ax2.set_xticks([])
        ax2.set_yticks([])
        data = ax2.imshow(cur[cluster_labels.argsort()].reshape((len(cluster_labels), 1)), cmap='jet', aspect='auto')
    # ax = plt.subplot(111)
    ax1.set_xticks(np.arange(len(gene_expression_top_var_headers_rows)))
    # ax1.set_yticks([])# [i for i in list(np.arange(len(tested_gene_expression_headers_columns))) if i % 10 == 0])
    ax1.set_xticklabels(e2g_convertor(gene_expression_top_var_headers_rows))
    ax1.set_yticklabels([])# [x for i, x in enumerate(cluster_labels[cluster_labels.argsort()]) if i % 10 == 0])
    for label in ax1.xaxis.get_ticklabels():
        label.set_fontsize(5)
        label.set_rotation(90)
    for label in ax1.yaxis.get_ticklabels():
        label.set_fontsize(3)
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    data = ax1.imshow(gene_expression_top_var[cluster_labels.argsort(), :], cmap='jet', aspect='auto')
    cb = plt.colorbar(data, ax=axes, fraction=0.05, pad=0.04)
    plt.savefig(os.path.join(constants.BASE_PROFILE, "output",
                             "heatmap_cluster_by_p_{}_{}_k={}_label_i={}_{}.png".format(constants.CANCER_TYPE,
                                                                          tested_gene_list_file_name.split(".")[0],
                                                                          n_clusters, label_index, time.time())))
    # "heatmap_cluster_by_p_{}_{}_k={}.svg".format(constants.CANCER_TYPE, tested_gene_list_file_name.split(".")[0], n_clusters)), format='svg', dpi=1200)
    # plt.clim(np.min(gene_expression_top_var),np.max(gene_expression_top_var))
    # plt.show()
    cb.remove()
    plt.cla()


def km_curve(labels_ids, survival_dataset, tested_gene_expression_headers_columns, gene_group , k=None, label_index=None):
    ax = plt.subplot(111)
    # box = ax.get_position()
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0 + box.height * 0.1,
    #                  box.width, box.height * 0.9])
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
    #       fancybox=True, shadow=True, ncol=5)

    kmf = KaplanMeierFitter()
    all_labels = np.array([y for x in labels_ids for y in x])
    for i, cur_labels in enumerate(labels_ids):
        # label_data = [[int(cur[4]), int(cur[5])] for cur in survival_dataset if cur[0] in cur_labels]
        # label_event, label_duration = zip(*label_data)
        # label_c_data = [[int(cur[4]), int(cur[5])] for cur in survival_dataset if cur[0] not in cur_labels]
        # label_event_c, label_duration_c = zip(*label_c_data)
        label_event = survival_dataset[np.in1d(survival_dataset[:, 0], cur_labels) & np.in1d(survival_dataset[:, 0], tested_gene_expression_headers_columns), 4].astype(np.int32)
        label_duration = survival_dataset[np.in1d(survival_dataset[:, 0], cur_labels) & np.in1d(survival_dataset[:, 0], tested_gene_expression_headers_columns), 3].astype(np.int32)

        # label_event= np.append(label_event, [0])
        # label_duration = np.append(label_duration, [543])

        labels_c = all_labels[~np.in1d(all_labels,cur_labels) & np.in1d(all_labels, tested_gene_expression_headers_columns)]
        label_event_c = survival_dataset[np.in1d(survival_dataset[:, 0], labels_c), 4].astype(np.int32)
        label_duration_c = survival_dataset[np.in1d(survival_dataset[:, 0], labels_c), 3].astype(np.int32)

        lr_results = logrank_test(label_duration, label_duration_c, label_event, label_event_c, alpha=.95)
        if len(label_duration) != 0:
            kmf.fit(list(label_duration), event_observed=list(label_event), label="cluster {} n={}, logrank pval = {}".format(i,len(label_duration), '{0:1.3e}'.format(lr_results.p_value))) # '%.7f' %
            kmf.plot(ax=ax)
            print "lrank: {}".format(lr_results.p_value)

    plt.ylim(0, 1);

    plt.title("clustering survival analysis");
    # plt.show()
    plt.savefig(os.path.join(constants.BASE_PROFILE,"output" ,"cluster_by_p_{}_{}_k={}_label_i={}_{}.png".format(constants.CANCER_TYPE, gene_group,k,label_index , time.time())))
    plt.cla()

def separate_headers(total_gene_expression):
    total_gene_expression = np.array(total_gene_expression)
    total_gene_expression_headers_columns = total_gene_expression[0][1:]
    total_gene_expression_headers_rows = total_gene_expression[1:, 0]
    total_gene_expression = total_gene_expression[1:]
    total_gene_expression = total_gene_expression[:, 1:]
    total_gene_expression = total_gene_expression.astype(np.float32)
    return total_gene_expression_headers_rows, total_gene_expression_headers_columns, total_gene_expression


def check_enrichment(gene_list):
    ensembl_for_url = re.sub("\.\d{1,2},", ",", gene_list)
    url = "http://david.abcc.ncifcrf.gov/api.jsp?type=ENSEMBL_GENE_ID&ids={}&tool=chartReport&annot=GOTERM_BP_DIRECT,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT,KEGG_PATHWAY".format(ensembl_for_url)
    return url

def print_to_excel(output_rows,gene_list_file_name,gene_expression_file_name,var_th_index):
    wb = Workbook()#ffff00
    ws = wb.active
    yellowFill = PatternFill(start_color='00FFFF00',
                          end_color='00FFFF00',
                          fill_type='solid')
    bd_regular = Side(style='thin', color="000000")
    border_regular = Border(left=bd_regular, top=bd_regular, right=bd_regular, bottom=bd_regular)

    bd_bold = Side(style='thick', color="000000")
    border_bold = Border(left=bd_bold, top=bd_bold, right=bd_bold, bottom=bd_bold)

    blueDarkFill = PatternFill(start_color='006699FF',
                             end_color='006699FF',
                             fill_type='solid')
    blueMediumFill = PatternFill(start_color='0099CCFF',
                               end_color='0099CCFF',
                               fill_type='solid')
    blueLightFill = PatternFill(start_color='00E6F3FF',
                                 end_color='00E6F3FF',
                                 fill_type='solid')
    border_regular = Border(left=bd_regular, top=bd_regular, right=bd_regular, bottom=bd_regular)


    headers = ["Description", "Genes", "URL", "GO_NS", "GO_ID", "GO_name", "nominal_pval", "FDR"]
    for i, header in enumerate(headers):
        ws['{}1'.format(chr(65+i))].border = border_regular
        ws['{}1'.format(chr(65+i))].fill = yellowFill
        ws['{}1'.format(chr(65+i))] = header
        ws.column_dimensions['{}'.format(chr(65 + i))].width = 30

    for k, cur in enumerate(headers):
        for i, cur in enumerate(output_rows):
            ws['{}{}'.format(chr(65+k), i+2)].border = border_regular
            ws['{}{}'.format(chr(65+k), i+2)].fill = blueLightFill
            ws['{}{}'.format(chr(65+k), i+2)] = cur[k]
            ws['{}{}'.format(chr(65+k), i+2)].alignment = Alignment(wrap_text=True)

    ws.column_dimensions["{}".format(chr(66+k))].width = 30
    ws["{}1".format(chr(66+k))].border = border_bold
    ws["{}1".format(chr(66+k))].fill = blueDarkFill
    ws["{}1".format(chr(66+k))] = "top var {}".format(var_th_index)
    wb.save(os.path.join(constants.OUTPUT_DIR,"CULSTERS_ENRICHMENT-{}-{}-topvar-{}-{}.xlsx".format(gene_list_file_name, gene_expression_file_name, var_th_index, time.time())))
