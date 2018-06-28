import re
import os
import math
from numpy import *
from decimal import Decimal, ROUND_HALF_EVEN
from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction
import numpy.random
from sklearn.datasets import fetch_mldata
import sklearn.preprocessing
import numpy as np
# import matplotlib as mpl
#mpl.use('Agg')
import scipy.special
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
# import matplotlib.pyplot as plt
# from matplotlib import style
# style.use("ggplot")
from sklearn import svm
from sklearn import svm
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.model_selection import PredefinedSplit
from sklearn.metrics import accuracy_score
import scipy
from scipy.stats import hypergeom
# from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import time
# from matplotlib.ticker import FormatStrFormatter
import math
import logging
from ge_pca_by_samples import genes_pca_by_samples

import json
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
#from clusters_similarity import *
from significant_enrichment import *
from gene_sets_correlation import *
from gene_sets_overlapping import *
#from svm import *
import fre
from fre import *
from feature_selection import *
from utils.param_builder import *
from genes_clustering import *
from patients_clustering import *
from cox_penotype_analysis import *
from cox_gene_analysis import *
from genes_statistic import *
from ge_prediction_by_mutations import *
from gene_rank_correlation import gene_correlation_scores
from genes_pca import genes_pca
from mutations_pca import mutation_pca
from mutations_pca_by_samples import mutation_pca_by_samples
from patients_clustring_filtered_by_mutations import cluster_patients_filtered_by_mutation
from mutations_pca_by_samples import plot_pca_by_samples
# RFE

## randomized by size
#
#  dataset = "SKCM"
# data_normalizaton = "counts"
# gene_expression_file_name, phenotype_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#
# gene_list_file_name = "protein_coding.txt"
# rounds=100
# rank_method = DISTANCE
# recursion_number_of_steps=1
# permutation=RANDOMIZED
# gene_sets_sizes = [13, 6, 7, 13]
# results = []
# for i in range(5):
#     for cur_size in gene_sets_sizes:
#         results.append(RFE(gene_list_file_name, gene_expression_file_name, phenotype_file_name, rank_method = rank_method, rounds=rounds, recursion_step_size=cur_size, recursion_number_of_steps=recursion_number_of_steps, pval_preprocessing_file_name = pval_preprocessing_file_name, permutation=permutation))
#     fre.print_to_excel(results=results, gene_sets_sizes=gene_sets_sizes,rank_method = rank_method, permutation=permutation)



##################### FEATURE SELECTION ################

# feature_selection(["rfe1000.txt"], "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv", gene_filter_file_name="protein_coding.txt", rounds=30, target_genes_subset = "mito.txt", recursion_step_size=1, feature_selection_method="rfe") #Primary Tumor , batches=100, epochs=20

##################### SETS CORRELATION ################

dataset = "UVM"
constants.update_dirs(CANCER_TYPE_u=dataset)
data_normalizaton = "fpkm_normalized_by_genes_standardization"
gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
find_sets_correlations(tested_gene_list_file_name=["mito.txt", "glycolysis_go.txt", "warburg_high.txt", "warburg_low.txt", "ldha_singular.txt"],
                       total_gene_list_file_name="protein_coding.txt", gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name)


##################### GENE CORRELATION SCORE ################

# gene_correlation_scores(tested_gene_list_file_name="mito.txt",
#                        total_gene_list_file_name="protein_coding.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", top_n=20)

##################### OVERLAPING ################

# find_sets_overlaps(tested_gene_list_file_name=["tca_pathcards.txt", "oxidative_phosphorylation_pathcards.txt", "glycolysis_pathcards.txt", "mito.txt", "oxidative.txt", "review.txt", "review_high.txt", "review_low.txt", "warburg_high.txt", "warburg_low.txt", "warburg.txt", "oxidative_HIF.txt", "pyruvate.txt", "ldha_singular.txt"],
#                        total_gene_list_file_name="protein_coding.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", phenotype_file_name="TCGA-SKCM.GDC_phenotype.tsv")


##################### CLUSTERS AND ENRICHMENT ################

# find_clusters_and_go_enrichment(tested_gene_list_file_name="rfe1000.txt", total_gene_list_file_name="protein_coding_long.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", phenotype_file_name="TCGA-SKCM.GDC_phenotype.tsv", var_th_index=None, start_k=2, end_k=6)

##################### CLUSTERS AND SURVIVAL ################


# # for x in range(10):
# #     random_set_file_name = generate_random_set(random_size=13, meta_gene_set="protein_coding_long.txt")
#
# for dataset in ["ESCA",  "ACC", "CHOL", "BLCA", "BRCA", "CESC", "COAD", "UCEC", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "DLBC", "LIHC", "LGG", "LUAD", "LUSC", "SKCM", "MESO", "UVM", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "STAD", "TGCT", "THYM", "THCA", "UCS"]:
# for dataset in ["SKCM"]:
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     tested_gene_list_file_name="mito.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=True
#     start_k=3
#     end_k=3
#     #[{"gender.demographic": {"type": "string", "value": ["male"]}},
#                                        # {"gender.demographic": {"type": "string", "value": ["female"]}}]
#     filter_expression = None
#     meta_groups = []
#     meta_groups=[
#         [
#             {"sample_type.samples": {"type": "string", "value": ["Primary Tumor", "Metastatic"]
#                                      },
#              "tumor_stage.diagnoses": {"type": "string", "value": ["stage 0", "stage i", "stage ia", "stage ib", "stage ic"]
#                                      }
#              },
#             {"sample_type.samples": {"type": "string", "value": ["Primary Tumor", "Metastatic"]
#                                      },
#              "tumor_stage.diagnoses": {"type": "string",
#                                        "value": ["stage iv", "stage iii", "stage iiia", "stage iiib", "stage iiic"]
#                                        }
#              }
#             ]
#     ]
#
#     # , "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]}
#
#     # for cur_tt in ["Primary Tumor"]:
#     filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor", "Metastatic"]},
#                           "tumor_stage.diagnoses": {"type": "string",
#                                                     "value": ["stage 0", "stage i", "stage ia", "stage ib", "stage ic", "stage iv", "stage iii", "stage iiia", "stage iiib",
#                                                               "stage iiic"]
#                                                     }
#                           }]
#     print "process {}".format(dataset)
#     find_clusters_and_survival(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups)

#################### SVM PREDICTION ##########################

# dataset = "SKCM"
# data_normalizaton = "counts"
# gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# # gene_list_file_names = ["warburg.txt", "warburg_high.txt", "warburg_low.txt", "mito_warburg.txt", "mito.txt", "test.txt", "test_2.txt", "review_high.txt", "review_low.txt", "review_low.txt", "review.txt"]
# gene_list_file_names = [ "kegg_melanoma.txt" ]
# gene_filter_file_name = "protein_coding.txt"
# rounds=100
# rank_method = DISTANCE
# results = []
# prediction_by_gene_expression(gene_list_file_names=gene_list_file_names, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, gene_filter_file_name=gene_filter_file_name, rounds=rounds, rank_method=rank_method, labels_permutation=constants.LABELS_NORMAL, compare_to_random = True,
#                               groups=[{"gender.demographic" : {"type": "string", "value" : ["male"]}},
#                                       {"gender.demographic": {"type": "string", "value": ["female"]}}])
#
# groups=[{"person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]}},
#                                       {"person_neoplasm_cancer_status": {"type": "string", "value": ["TUMOR FREE"]}}]

#################### COX PHENOTYPE ANALYSIS ##########################

# dataset = "SKCM"
# data_normalizaton = "counts"
# gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# # gene_list_file_names = ["warburg.txt", "warburg_high.txt", "warburg_low.txt", "mito_warburg.txt", "mito.txt", "test.txt", "test_2.txt", "review_high.txt", "review_low.txt", "review_low.txt", "review.txt"]
# test_independently = False
# filtered_out = ['bcr', 'days_to_initial_pathologic_diagnosis', 'disease_code', 'informed_consent_verified', 'project_code', 'withdrawn', 'year_of_dcc_upload', "therapy_type_notes", "pathology_report_file_name", "tx_on_clinical_trial", "units", "year_of_form_completion", "vial_number", "ethnicity.demographic", "is_ffpe.samples", "project.tissue_source_site", "name.tissue_source_site","bcr_id.tissue_source_site", "project_id.project", "submitter_id", "therapy_type", "therapy_type_notes", 'classification_of_tumor.diagnoses', 'last_known_disease_status.diagnoses', 'prior_malignancy.diagnoses', 'progression_or_recurrence.diagnoses', 'tumor_grade.diagnoses', 'disease_type', 'primary_site', 'disease_type.project', 'name.project', 'primary_site.project', "vital_status.diagnoses", "day_of_dcc_upload", "month_of_dcc_upload", "tissue_retrospective_collection_indicator", "tissue_prospective_collection_indicator", "code.tissue_source_site", "tissue_source_site", "site_of_resection_or_biopsy.diagnoses", "disease_type", "primary_site", "disease_type.project", "name.project", "primary_site.project", "project_id.project", "submitter_id", "bcr_id.tissue_source_site", "code.tissue_source_site", "name.tissue_source_site", "project.tissue_source_site", "is_ffpe.samples", "oct_embedded.samples", "sample_type_id.samples", "state.samples", "system_version",  "age_at_diagnosis.diagnoses", "day_of_form_completion", "patient_id", "batch_number"]
# filtered_in = ["sample_type.samples", 'gender.demographic', "pathologic_M","pathologic_T","pathologic_N", "tumor_stage.diagnoses", "melanoma_clark_level_value", "days_to_birth.diagnoses" ,"weight.exposures","height.exposures", "person_neoplasm_cancer_status"]
# filter_type = FILTER_IN
# filter_na_by_rows = True
# f = file(os.path.join(constants.OUTPUT_DIR,"output_expression"), 'w+')
# gene_expression_file_name, phenotype_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# cox_phenotype(test_independently=test_independently, filtered_out=filtered_out, filtered_in=filtered_in, filter_type=filter_type, filter_na_by_rows=filter_na_by_rows, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name)
# f.close()

#################### COX GENE ANALYSIS ##########################
# dataset = "ACC"
# data_normalizaton = "counts"
# gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
# constants.update_dirs(CANCER_TYPE_u=dataset)
# test_independently = False
# filtered_out = []
# filtered_in = []
# tested_genes_file_name = "mito.txt"
# filter_type = FILTER_IN
# filter_na_by_rows = True
# filter_expression = None
# f = file(os.path.join(constants.BASE_PROFILE, "output","cox_genes_{}_{}.txt".format(dataset, time.time())), 'w+')
# # filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor"]}, "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR", "TUMOR FREE", ""]}}]
# results = cox_gene(test_independently, filtered_out, filtered_in, filter_type, filter_na_by_rows,tested_genes_file_name = tested_genes_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, filter_expression=filter_expression)
# f.write(results)
# f.close()


##################### CLUSTERS AND SURVIVAL ################

# for dataset in ["PAAD"]:
#     meta_groups=[]
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     tested_gene_list_file_name="mito.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=True
#     start_k=3
#     end_k=3
#     filter_expression = None
#     d_codes = []
#     # for cur in meta_groups[0]:
#     #     d_codes.append(cur["disease_code"]["value"][0])
#
#     # for cur_tt in ["Primary Tumor"]:
#     filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor","Metastatic"]},
#                           # "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]} ,
#                           # "disease_code": {"type": "string",
#                           #              "value": d_codes
#                           #              }
#                           }]
#     print "process {}".format(dataset)
#     find_clusters_and_survival(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups)

############# PCA for gene expression ##################

# for dataset in ["PANCAN"]:
#     meta_groups= [] # json.load("utils/tumor_type_groups.json")
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     tested_gene_list_file_name="mito.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=True
#     start_k=2
#     end_k=2
#     #[{"gender.demographic": {"type": "string", "value": ["male"]}},
#                                        # {"gender.demographic": {"type": "string", "value": ["female"]}}]
#     filter_expression = None
#
#     d_codes = []
#     for cur in meta_groups[0]:
#         d_codes.append(cur["disease_code"]["value"][0])
#
#     # , "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]}
#
#     # for cur_tt in ["Primary Tumor"]:
#     filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor","Metastatic"]},
#                           # "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]} ,
#                           # "disease_code": {"type": "string",
#                           #              "value": d_codes
#                           #              }
#                           }]
#     print "process {}".format(dataset)
#     genes_pca(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups)


############ PREDICTION_BY_MUTATION ###################

# for dataset in ["UVM"]:
#     if dataset == "PANCAN": continue
#     meta_groups = None
#     meta_groups= [json.load(file("groups/mutation_binary.json"))]
#
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "fpkm"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     random_set_file_name = generate_random_set(random_size=25, meta_gene_set="protein_coding_long.txt")
#     tested_gene_list_file_name= "glycolysis_go.txt" # random_set_file_name #
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=False
#     integ=False
#     start_k=2
#     end_k=2
#     omitted_genes = [] # ["TTN", "BRAF", "GNAQ", "GNA11"]
#     filter_expression = None
#     d_codes = []
#     # if meta_groups is not None:
#     #     for cur in meta_groups[0]:
#     #         d_codes.append(cur["tumor_stage.diagnoses"]["value"][0])
#     # , "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]}
#
#     # for cur_tt in ["Primary Tumor"]:
#     filter_expression = [json.load(file("filters/mutation_binary.json"))]
#     print "process {}".format(dataset)
#     phenotype_labels_heatmap = None#["breslow_depth_value"]
#     # find_clusters_and_survival(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups)
#     predict_ge_by_mutation(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups, phenotype_labels_heatmap=phenotype_labels_heatmap, integ=integ, min_ratio=0.4, omitted_genes=omitted_genes)
#     #cluster_patients_filtered_by_mutation(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups, integ=True ,min_ratio=0.4)

############ mutation PCA ##################

# for dataset in ["PANCAN"]:
#     meta_groups = json.load("utils/tumor_type_groups.json")
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     tested_gene_list_file_name="mito.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=True
#     start_k=2
#     end_k=2
#     filter_expression = None
#     d_codes = []
#     # for cur in meta_groups[0]:
#     #     d_codes.append(cur["disease_code"]["value"][0])
#
#     filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor","Metastatic"]},
#                           # "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]} ,
#                           # "disease_code": {"type": "string",
#                           #              "value": d_codes
#                           #              }
#                           }]
#     print "process {}".format(dataset)
#     mutation_pca(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups)

############ mutation PCA ##################

# for dataset in ["UVM"]:
#     # meta_groups = json.load("utils/tumor_type_groups.json")
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     tested_gene_list_file_name="mito.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=True
#     start_k=2
#     end_k=2
#     filter_expression = None
#     d_codes = []
#     # for cur in meta_groups[0]:
#     #     d_codes.append(cur["disease_code"]["value"][0])
#
#     filter_expression = [{"sample_type.samples": {"type": "string", "value": ["Primary Tumor"]},
#                           # "person_neoplasm_cancer_status" : {"type": "string", "value" : ["WITH TUMOR"]} ,
#                           # "disease_code": {"type": "string",
#                           #              "value": d_codes
#                           #              }
#                           }]
#     print "process {}".format(dataset)
#     mutation_pca_by_samples(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = None)

############ PREDICTION_BY_MUTATION ###################

# for dataset in ["UVM"]:
#     if dataset == "PANCAN": continue
#     meta_groups = None
#     meta_groups= [json.load(file("groups/mutation_cluster.json"))]
#
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     random_set_file_name = generate_random_set(random_size=25, meta_gene_set="protein_coding_long.txt")
#     tested_gene_list_file_name= "warburg_high.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     filter_expression = None
#     print "process {}".format(dataset)
#     genes_pca_by_samples(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, filter_expression= filter_expression, meta_groups = meta_groups)
#     # predict_ge_by_mutation(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups, phenotype_labels_heatmap=phenotype_labels_heatmap, integ=integ, min_ratio=0.9, omitted_genes=omitted_genes)
#     # cluster_patients_filtered_by_mutation(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, mutation_file_name=mutation_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, filter_expression= filter_expression, meta_groups = meta_groups, integ=True ,min_ratio=0.4)
