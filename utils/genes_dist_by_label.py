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
    data_primary = None
    data_metastatic = None
    labels = []
    primary_labeled, metastatic_labeled = load_expression_profile_by_gene_and_tumor_type(tested_gene_file_name, expression_profile_file_name, phenotype_file_name, gene_filter_file_name)

    for cur in primary_labeled[1:]:
        if data_primary is None:
            data_primary = np.array(cur[1:])
        else:
            data_primary = np.vstack((data_primary,cur[1:]))
        #data_primary.append(cur[1:].astype(np.float))
    data_primary = np.array([cur_j for cur_i in data_primary for cur_j in cur_i]).astype(np.float)
    # print scipy.stats.mstats.normaltest(data_primary)
    sns.distplot(data_primary)
    # plt.savefig("primary_dist_{}".format(expression_profile_file_name.split("_")[1].split('.')[0]))

    for cur in metastatic_labeled[1:]:
        if data_metastatic is None:
            data_metastatic = np.array(cur[1:])
        else:
            data_metastatic = np.vstack((data_metastatic, cur[1:]))
        # data_metastatic.append(cur[1:].astype(np.float))
    data_metastatic = np.array([cur_j for cur_i in data_metastatic for cur_j in cur_i]).astype(np.float)
    # print scipy.stats.mstats.normaltest(data_metastatic)
    sns.distplot(data_metastatic)
    plt.savefig(os.path.join(OUTPUT_DIR, "dist_comparison_{}".format(expression_profile_file_name.split("_")[1].split('.')[0])))
    plt.cla()

load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_fpkm.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_fpkm-uq.tsv", "TCGA-SKCM.GDC_phenotype.tsv")

