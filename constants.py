USE_CACHE = True
PHENOTYPE_FORMAT = "GDC"
BASE_PROFILE="D:\\omics\\"
BASE_DATASET= "{}GDC-TCGA\\melanoma\\".format(BASE_PROFILE)
CACHE_DIR = "{}cache\\".format(BASE_DATASET)
DICT_DIR = "{}dictionaries\\".format(BASE_DATASET)
OUTPUT_DIR = "{}output\\".format(BASE_DATASET)
LIST_DIR = "{}list\\".format(BASE_DATASET)
TCGA_DATA_DIR = "{}tcga_data\\".format(BASE_DATASET)

SEPARATOR = "@%@"

LABEL_ID = "sample_type.samples"
PRIMARY_TUMOR = "Primary Tumor"
METASTATIC = "Metastatic"

LABELS_NORMAL = "labels_normal"
LABELS_SHUFFLE = "labels_shuffle"
LABELS_RANDOM = "labels_random"
LABELS_ALTERNATED = "labels_alternated"
LABELS_INVERTED = "labels_inverted"

def update_dirs(BASE_DIR="D:\\omics\\", DATASET_DIR="GDC-TCGA\\melanoma\\"):
    global BASE_PROFILE
    global CACHE_DIR
    global OUTPUT_DIR
    global LIST_DIR
    global TCGA_DATA_DIR
    global DICTIONARIES_DIR
    BASE_PROFILE=BASE_DIR
    BASE_DATASET= "{}{}".format(BASE_PROFILE,DATASET_DIR)
    CACHE_DIR= "{}cache\\".format(BASE_DATASET)
    OUTPUT_DIR = "{}output\\".format(BASE_DATASET)
    LIST_DIR = "{}list\\".format(BASE_DATASET)
    TCGA_DATA_DIR = "{}tcga_data\\".format(BASE_DATASET)
    DICTIONARIES_DIR = "{}dictionaries\\".format(BASE_DATASET)

update_dirs()