import scipy.special
import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
import scipy
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from scipy.stats import pearsonr
from utils.ensembl2gene_symbol import e2g_convertor
############################ (2) correlations #############################


# () main
def gene_correlation_scores(tested_gene_list_file_name, total_gene_list_file_name, gene_expression_file_name, gene_filter_file_name=None, top_n=2000):
    print "about ot analyse: {}".format(tested_gene_list_file_name)
    # fetch gene expression by gene_id, divided by tumor type
    total_gene_expression = np.array(load_gene_expression_profile_by_genes(total_gene_list_file_name, gene_expression_file_name, gene_filter_file_name))
    gene_ids = load_gene_list(tested_gene_list_file_name)
    ranks_dict = {}
    ranks_score = []
    ranks = []
    for cur in gene_ids:
        ranks.append([])
        cur_expression = total_gene_expression[np.where(total_gene_expression[:,0]==cur)][0]
        for cur_prot_expression in total_gene_expression[1:]:
            prs =pearsonr(cur_prot_expression[1:].astype(np.float32), cur_expression[1:].astype(np.float32))[0]
            if not math.isnan(prs) and cur_prot_expression[0] not in gene_ids:
                ranks[-1].append((cur_prot_expression[0],abs(prs), prs>0))
        ranks[-1] = sorted(ranks[-1], key=lambda x: x[1], reverse=True)[:top_n]

        for cur in ranks[-1]:
            if not ranks_dict.has_key(cur[0]):
                ranks_dict[cur[0]] = []
            ranks_dict[cur[0]].append(cur[1])

    for k,v in ranks_dict.iteritems():
        ranks_score.append((k,sum(v), e2g_convertor([k])[0]))
    ranks_score = sorted(ranks_score, key=lambda x: x[1], reverse=True)[:top_n]
    print ranks_score

    f = file(os.path.join(constants.LIST_DIR, "corr_{}_top_{}.txt".format(tested_gene_list_file_name.split(".")[0], top_n)), 'w+')
    f.write("\r\n".join([x[0] for x in ranks_score]))
    f.close()
