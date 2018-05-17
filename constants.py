USE_CACHE = True
BASE_PROFILE="D:\\omics\\"
BASE_MELANOMA= "{}GDC-TCGA\\melanoma\\".format(BASE_PROFILE)
CACHE_DIR = "{}cache\\".format(BASE_MELANOMA)
OUTPUT_DIR = "{}output\\".format(BASE_MELANOMA)
LIST_DIR = "{}list\\".format(BASE_MELANOMA)
TCGA_DATA_DIR = "{}tcga_data\\".format(BASE_MELANOMA)
SEPARATOR = "@%@"

LABEL_ID = "sample_type.samples"
PRIMARY_TUMOR = "Primary Tumor"
METASTATIC = "Metastatic"
