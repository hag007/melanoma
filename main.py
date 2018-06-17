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

##################### CORRELATION ################

# find_sets_correlations(tested_gene_list_file_name=["tca_pathcards.txt", "oxidative_phosphorylation_pathcards.txt", "glycolysis_pathcards.txt", "mito.txt", "oxidative.txt", "review.txt", "review_high.txt", "review_low.txt", "warburg_high.txt", "warburg_low.txt", "warburg.txt", "oxidative_HIF.txt", "pyruvate.txt", "ldha_singular.txt"],
#                        total_gene_list_file_name="protein_coding.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", phenotype_file_name="TCGA-SKCM.GDC_phenotype.tsv")

##################### OVERLAPING ################

find_sets_overlaps(tested_gene_list_file_name=["tca_pathcards.txt", "oxidative_phosphorylation_pathcards.txt", "glycolysis_pathcards.txt", "mito.txt", "oxidative.txt", "review.txt", "review_high.txt", "review_low.txt", "warburg_high.txt", "warburg_low.txt", "warburg.txt", "oxidative_HIF.txt", "pyruvate.txt", "ldha_singular.txt"],
                       total_gene_list_file_name="protein_coding.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", phenotype_file_name="TCGA-SKCM.GDC_phenotype.tsv")


##################### CLUSTERS AND ENRICHMENT ################

# find_clusters_and_go_enrichment(tested_gene_list_file_name="rfe1000.txt", total_gene_list_file_name="protein_coding_long.txt", gene_expression_file_name="TCGA-SKCM.htseq_counts.tsv", phenotype_file_name="TCGA-SKCM.GDC_phenotype.tsv", var_th_index=None, start_k=2, end_k=6)

##################### CLUSTERS AND SURVIVAL ################

# for dataset in ["ACC","PAAD", "LAML", "CHOL", "UVM", "LUAD", "LUSC", "BRCA"]:
# for dataset in ["SKCM"]:
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "counts_normalized_by_patients"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, pval_preprocessing_file_name = build_gdc_params(dataset=dataset, data_normalizaton=data_normalizaton)
#     tested_gene_list_file_name="mito.txt"
#     total_gene_list_file_name="protein_coding_long.txt"
#     var_th_index=None
#     is_unsupervised=True
#     start_k=2
#     end_k=2
#     groups=None #[{"gender.demographic": {"type": "string", "value": ["male"]}},
#                                        # {"gender.demographic": {"type": "string", "value": ["female"]}}]
#
#     find_clusters_and_survival(tested_gene_list_file_name=tested_gene_list_file_name, total_gene_list_file_name=total_gene_list_file_name, gene_expression_file_name=gene_expression_file_name, phenotype_file_name=phenotype_file_name, survival_file_name=survival_file_name, var_th_index=var_th_index, is_unsupervised=is_unsupervised, start_k=start_k, end_k=end_k, groups=groups)

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

# test_independently = False
# filtered_out = []
# filtered_in = ["ENSG00000004975.10", "ENSG00000010404.16", "ENSG00000016391.9", "ENSG00000016490.14", "ENSG00000042980.11", "ENSG00000049249.7", "ENSG00000052126.13", "ENSG00000057019.14", "ENSG00000064115.9", "ENSG00000065615.12", "ENSG00000068985.4", "ENSG00000074201.7", "ENSG00000075651.14", "ENSG00000083290.18", "ENSG00000084093.14", "ENSG00000086544.2", "ENSG00000087116.12", "ENSG00000088367.19", "ENSG00000088682.12", "ENSG00000100077.13", "ENSG00000100167.18", "ENSG00000101310.13", "ENSG00000104047.13", "ENSG00000104133.13", "ENSG00000105137.11", "ENSG00000105383.13", "ENSG00000106538.8", "ENSG00000108242.11", "ENSG00000109099.12", "ENSG00000109181.10", "ENSG00000113407.12", "ENSG00000115216.12", "ENSG00000118518.14", "ENSG00000120458.8", "ENSG00000121741.15", "ENSG00000124098.9", "ENSG00000124104.17", "ENSG00000124140.11", "ENSG00000124422.10", "ENSG00000124875.8", "ENSG00000125779.20", "ENSG00000125995.14", "ENSG00000126231.12", "ENSG00000129128.11", "ENSG00000129353.13", "ENSG00000129450.7", "ENSG00000129451.10", "ENSG00000131738.8", "ENSG00000132155.10", "ENSG00000132694.17", "ENSG00000133475.15", "ENSG00000134755.13", "ENSG00000134759.12", "ENSG00000134760.5", "ENSG00000134762.15", "ENSG00000136153.18", "ENSG00000136682.13", "ENSG00000136694.8", "ENSG00000137434.10", "ENSG00000137642.11", "ENSG00000137693.12", "ENSG00000137707.12", "ENSG00000137968.15", "ENSG00000138246.14", "ENSG00000140199.10", "ENSG00000140274.12", "ENSG00000141577.12", "ENSG00000142252.9", "ENSG00000142621.18", "ENSG00000143515.15", "ENSG00000143520.6", "ENSG00000143554.12", "ENSG00000144451.17", "ENSG00000145920.13", "ENSG00000147669.9", "ENSG00000147687.15", "ENSG00000148344.10", "ENSG00000153292.14", "ENSG00000153790.10", "ENSG00000154222.13", "ENSG00000156453.12", "ENSG00000158050.4", "ENSG00000158122.10", "ENSG00000158773.13", "ENSG00000159335.14", "ENSG00000159450.11", "ENSG00000159496.13", "ENSG00000159516.8", "ENSG00000161791.12", "ENSG00000162068.1", "ENSG00000162365.10", "ENSG00000163191.5", "ENSG00000163206.5", "ENSG00000163207.6", "ENSG00000163214.19", "ENSG00000163217.1", "ENSG00000163431.12", "ENSG00000164366.3", "ENSG00000164694.15", "ENSG00000164823.8", "ENSG00000164853.8", "ENSG00000165792.16", "ENSG00000165795.19", "ENSG00000165949.11", "ENSG00000166181.11", "ENSG00000166394.13", "ENSG00000166532.14", "ENSG00000166669.12", "ENSG00000166866.11", "ENSG00000167642.11", "ENSG00000167654.16", "ENSG00000167740.8", "ENSG00000167751.11", "ENSG00000167754.11", "ENSG00000167755.12", "ENSG00000167757.12", "ENSG00000167767.12", "ENSG00000167768.4", "ENSG00000168140.4", "ENSG00000168702.15", "ENSG00000168906.11", "ENSG00000169032.8", "ENSG00000169446.5", "ENSG00000169469.8", "ENSG00000169508.6", "ENSG00000169548.3", "ENSG00000169592.13", "ENSG00000170421.10", "ENSG00000170425.3", "ENSG00000170448.10", "ENSG00000170464.8", "ENSG00000170476.14", "ENSG00000170782.3", "ENSG00000171121.15", "ENSG00000171345.12", "ENSG00000171396.11", "ENSG00000171402.13", "ENSG00000171953.14", "ENSG00000172476.3", "ENSG00000172548.13", "ENSG00000172817.3", "ENSG00000172845.12", "ENSG00000174948.5", "ENSG00000175115.10", "ENSG00000175311.6", "ENSG00000175701.9", "ENSG00000175792.10", "ENSG00000176087.13", "ENSG00000176182.5", "ENSG00000178171.9", "ENSG00000178358.4", "ENSG00000178363.4", "ENSG00000178591.6", "ENSG00000178917.13", "ENSG00000178928.7", "ENSG00000179144.4", "ENSG00000179172.8", "ENSG00000179476.6", "ENSG00000181323.7", "ENSG00000182035.10", "ENSG00000182584.4", "ENSG00000183346.6", "ENSG00000183753.8", "ENSG00000184144.8", "ENSG00000184321.1", "ENSG00000184454.6", "ENSG00000184730.9", "ENSG00000185477.4", "ENSG00000185634.10", "ENSG00000185963.12", "ENSG00000186472.18", "ENSG00000186827.9", "ENSG00000186838.12", "ENSG00000186844.5", "ENSG00000187170.4", "ENSG00000187175.5", "ENSG00000187210.11", "ENSG00000188060.6", "ENSG00000188086.11", "ENSG00000188095.4", "ENSG00000188290.9", "ENSG00000188372.13", "ENSG00000188505.4", "ENSG00000188909.4", "ENSG00000188997.6", "ENSG00000189181.4", "ENSG00000189269.11", "ENSG00000189326.4", "ENSG00000189430.11", "ENSG00000196539.3", "ENSG00000196730.11", "ENSG00000196800.5", "ENSG00000197081.11", "ENSG00000197635.8", "ENSG00000198064.12", "ENSG00000198853.10", "ENSG00000203780.9", "ENSG00000203782.5", "ENSG00000203784.2", "ENSG00000203785.7", "ENSG00000203818.6", "ENSG00000204420.7", "ENSG00000204536.12", "ENSG00000204538.3", "ENSG00000204859.10", "ENSG00000205060.9", "ENSG00000205413.6", "ENSG00000206072.11", "ENSG00000206073.9", "ENSG00000215845.9", "ENSG00000223569.4", "ENSG00000235863.3", "ENSG00000237254.2", "ENSG00000240382.3", "ENSG00000241755.1", "ENSG00000244067.2", "ENSG00000256806.3", "ENSG00000261147.1", "ENSG00000264187.1", "ENSG00000265118.4", "ENSG00000268089.2", "ENSG00000269720.1", "ENSG00000281899.1"]
# filter_type = FILTER_IN
# filter_na_by_rows = True
# f = file(os.path.join(constants.OUTPUT_DIR,"output_expression"), 'w+')
# results = cox_gene(test_independently, filtered_out, filtered_in, filter_type, filter_na_by_rows)
# f.close()


