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
            if cur_group is None:
                groups[i] = np.array(cur[1:])
            else:
                groups[i] = np.vstack((groups[i],cur[1:]))
        data_primary = np.array([cur_j for cur_i in groups[i] for cur_j in cur_i]).astype(np.float)
        sns.distplot(data_primary)
        plt.savefig(os.path.join(OUTPUT_DIR, "dist_comparison_{}".format(expression_profile_file_name.split("_")[1].split('.')[0])))
        plt.cla()

load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_counts.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_fpkm.tsv", "TCGA-SKCM.GDC_phenotype.tsv")
load_gene_expression_by_label("protein_coding.txt", "TCGA-SKCM.htseq_fpkm-uq.tsv", "TCGA-SKCM.GDC_phenotype.tsv")

