from matplotlib import style
style.use("ggplot")
import seaborn as sns
sns.set(color_codes=True)
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *

DISTANCE = "distance"
LOGISTIC_REGRESSION = "logistic_regression"

def load_gene_expression_by_label(tested_gene_file_name, expression_profile_file_name, phenotype_file_name, gene_filter_file_name=None):
    dist_data = []
    groups = load_expression_profile_by_labelling(tested_gene_file_name, expression_profile_file_name, phenotype_file_name, gene_filter_file_name)
    for cur in groups:
        dist_data.append(None)
    for i, cur_group in enumerate(groups):
        for cur in cur_group[1:]:
            if dist_data[i] is None:
                dist_data[i] = np.array(cur[1:])
            else:
                dist_data[i] = np.vstack((dist_data[i],cur[1:]))
        data_primary = np.array([cur_j for cur_i in dist_data[i] for cur_j in cur_i]).astype(np.float)
        sns.distplot(data_primary)
        plt.savefig(os.path.join(OUTPUT_DIR, "dist_comparison_{}_group_{}".format(expression_profile_file_name.split("_")[1].split('.')[0],i)))
        plt.cla()

load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_counts_normalized_by_patients_standardization.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_fpkm_normalized_by_patients_standardization.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_fpkm-uq_normalized_by_patients_standardization.tsv", "TCGA-SKCM.GDC_phenotype.tsv")

# load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_counts_normalization_.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
# load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_fpkm.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
# load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_fpkm-uq.tsv", "TCGA-SKCM.GDC_phenotype.tsv")



