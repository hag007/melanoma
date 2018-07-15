import gzip
import shutil
# import math
# import os
# import pandas as pd
# import patsy
# import sys
# import numpy.linalg as la
# import numpy as np
# import constants
# from utils.param_builder import build_gdc_params

with gzip.open("GDC-PANCAN.GDC_phenotype.tsv.gz", 'rb') as f_in:
    with open("GDC-PANCAN.GDC_phenotype.tsv", 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


# pheno = []
# dat  =[]
# for dataset in ["STAD", "PAAD", "UVM", "BRCA"]:
#
#     print "current dataset: {}".format(dataset)
#     constants.update_dirs(CANCER_TYPE_u=dataset)
#     data_normalizaton = "fpkm"
#     gene_expression_file_name, phenotype_file_name, survival_file_name, mutation_file_name, mirna_file_name, pval_preprocessing_file_name = build_gdc_params(
#         dataset=dataset, data_normalizaton=data_normalizaton)
#
#     # pheno.append(pd.read_table(os.path.join(constants.TCGA_DATA_DIR, phenotype_file_name), index_col=0))
#     dat.append(pd.read_table(os.path.join(constants.TCGA_DATA_DIR, gene_expression_file_name), index_col=0))
#
# # pd.concat(pheno, axis=0).to_csv("pheno_bc", sep="\t")
# dat = pd.concat(dat, axis=1)
# dat = dat.fillna(0)
# dat.index.name = 'ENSEMBL_ID'
# dat.to_csv("dat", sep="\t")