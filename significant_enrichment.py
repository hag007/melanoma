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

############################ (2) significance expression and proportion differntiations #############################


def plot_pvalues(y_axis, x_axis, th,output_file_name, is_bigger_better=False):
    n, bins, patches = mpl.pyplot.hist(y_axis, x_axis)
    for c, p in zip(bins, patches):
        if is_bigger_better:
            th_condition = c < th
        else:
            th_condition = c > th

        if th_condition:
            color = 'blue'
        else:
            color = 'red'
        plt.setp(p, 'facecolor', color)
    plt.savefig(os.path.join(OUTPUT_DIR, output_file_name))
    plt.cla()


def plot_pvalues_log_scaled(y_axis, x_axis, th,output_file_name):
    plot_pvalues([-math.log(cur, 10) for cur in y_axis], x_axis, -math.log(th, 10), output_file_name, is_bigger_better=True)


def summarize_genes_proportion(tested_gene_list, total_gene_list, gene_pval_pair, true_counter):
    significant_tested_gene_list = [((i + 1), cur[0], cur[1]) for i, cur in enumerate(gene_pval_pair) if
                                    cur[0] in tested_gene_list and i < true_counter]
    included_tested_genes = [cur for cur in tested_gene_list if
                             cur in total_gene_list or cur[:cur.find('.')] in total_gene_list]
    included_tested_genes_size = len(included_tested_genes)
    significant_tested_gene_list_size = len(significant_tested_gene_list)
    print "total tested genes in true hypothsis: {} out of possible {}".format(significant_tested_gene_list_size, included_tested_genes_size)
    tested_gene_list_size = len(tested_gene_list)
    total_gene_list_size = len(total_gene_list)
    results_table = []
    expected_actual_difference_list = []
    expected_actual_ratio_difference_list = []
    rank_in_total_list_list = []
    z_test_proportion_test_list = []
    for i, cur in enumerate(significant_tested_gene_list):
        rank_in_tesed_list = (i + 1)
        rank_in_total_list = cur[0]
        ensembel_id = cur[1]
        p_val = cur[2]
        expected_quantity = rank_in_total_list * (included_tested_genes_size / (total_gene_list_size * 1.0))
        expected_proportion = included_tested_genes_size / (total_gene_list_size * 1.0)
        actual_quantity = (i + 1)
        actual_proportion = (i + 1) / (cur[0] * 1.0)
        expected_actual_difference = actual_quantity - expected_quantity
        expected_actual_ratio_difference = expected_actual_difference / (expected_quantity * 1.0)
        z_test_proportion_test = (actual_proportion - expected_proportion) / math.sqrt(
            (expected_proportion * (1 - expected_proportion)) / rank_in_total_list)
        results_table.append([rank_in_tesed_list, rank_in_total_list, ensembel_id, p_val,
                              expected_quantity, expected_proportion,
                              actual_quantity, actual_proportion,
                              expected_actual_difference, expected_actual_ratio_difference,
                              z_test_proportion_test])

        expected_actual_difference_list.append(expected_actual_difference)
        expected_actual_ratio_difference_list.append(expected_actual_ratio_difference)
        z_test_proportion_test_list.append(z_test_proportion_test)
        rank_in_total_list_list.append(rank_in_total_list)
    return expected_actual_difference_list, expected_actual_ratio_difference_list, rank_in_total_list_list, z_test_proportion_test_list


def plot_genes_proportion(expected_actual_difference_list, expected_actual_ratio_difference_list, z_test_proportion_test_list, rank_in_total_list_list, total_significant_hypotheses_size, expected_tested_genes_ratio, tested_gene_list_file_name):
    z_score_threshold_two_way = 1.96
    tested_genes_size = len(rank_in_total_list_list)
    y_counter = [min(i, tested_genes_size) for i in range(1,tested_genes_size+1)]
    plt.plot(rank_in_total_list_list[1:], y_counter[1:], label="number of significant values (n)")
    plt.plot(rank_in_total_list_list[1:], expected_actual_difference_list[1:], label="actual-expected significant hypo. difference (n)")
    plt.plot([total_significant_hypotheses_size, total_significant_hypotheses_size], [-20, tested_genes_size + 5], label="True hypotheses threshold")
    plt.plot([0, total_significant_hypotheses_size], [0, tested_genes_size],
             color="gray")
    plt.plot([0, total_significant_hypotheses_size], [0, expected_tested_genes_ratio],
             color="black")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "{}_sum_n".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')])))
    plt.cla()
    plt.plot(rank_in_total_list_list[1:], expected_actual_ratio_difference_list[1:], label="actual/expected proportion ratio")
    plt.plot(rank_in_total_list_list[1:], z_test_proportion_test_list[1:], label="z_score")
    plt.plot(rank_in_total_list_list[1:], [z_score_threshold_two_way for i in range(1,tested_genes_size)], label="z_score threshold (two-way)")
    plt.plot([total_significant_hypotheses_size, total_significant_hypotheses_size], [-0.5, 3.5], label="True hypotheses threshold")
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "{}_sum_p".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')])))

# mHGT DP
def calc_num_of_non_extremer_paths(non_extremer_paths_DP_table, HGTs, mHGT, n, b):
    if n==0 and b==0:
        # non_extremer_paths_DP_table[b][n] = 1
        return 1
    elif b==-1 or b>n or (HGTs[b][n] < mHGT and (b<len(HGTs)-1 or n<len(HGTs[0])-1)):
        # non_extremer_paths_DP_table[b][n] = 0
        return 0
    elif non_extremer_paths_DP_table[b][n] == -1:
        non_extremer_paths_DP_table[b][n] = long(calc_num_of_non_extremer_paths(non_extremer_paths_DP_table, HGTs, mHGT, n-1, b)) + long(calc_num_of_non_extremer_paths(non_extremer_paths_DP_table, HGTs, mHGT, n-1, b-1))


    return non_extremer_paths_DP_table[b][n]

# (2) main
def find_expression_significance(tested_gene_list_file_name, total_gene_list_file_name, gene_expression_file_name, phenotype_file_name, gene_filter_file_name=None, tested_gene_list_path=None, total_gene_list_path=None, gene_expression_path=None, phenotype_path=None, gene_filter_file_path=None, source="GDC-TCGA",dataset="melanoma", hgt_preprocessing_file_name = None, pval_preprocessing_file_name= None, N = None, B = None):
    print "about ot analyse: {}".format(tested_gene_list_file_name)
    # fetch gene expression by gene_id, divided by tumor type
    groups = load_expression_profile_by_labelling(total_gene_list_file_name, gene_expression_file_name, phenotype_file_name, gene_filter_file_name, total_gene_list_path, gene_expression_path, phenotype_path, gene_filter_file_path, source, dataset)
    group_0_expression = groups[0]
    group_1_expression = groups[1]
    group_0_expression = np.rot90(np.flip(group_0_expression, 1), k=-1, axes=(1,0))
    group_1_expression = np.rot90(np.flip(group_1_expression, 1), k=-1, axes=(1, 0))

    # test pval for significance differentiation between label values (primar vs metastatic)
    if os.path.isfile(os.path.join(CACHE_DIR, pval_preprocessing_file_name)) and USE_CACHE:
        gene_pval_pair = load_sets(os.path.join(CACHE_DIR, pval_preprocessing_file_name))
        print "pval loaded from file"
    else:
        pvals = []
        gene_symbols = []
        for  i in range(1,len(group_0_expression)):
            cur_pval = scipy.stats.ttest_ind([float(c) for c in group_0_expression[i][1:]], [float(c) for c in group_1_expression[i][1: ]])[1]
            if not math.isnan(cur_pval):
                pvals.append(cur_pval)
                gene_symbols.append(group_0_expression[i][0])

        # sort gene_id-pval pairs by pval
        gene_pval_pair = zip(gene_symbols,pvals)
        gene_pval_pair.sort(key=lambda x: x[1], reverse=False)
        save_sets(gene_pval_pair, os.path.join(CACHE_DIR, os.path.join(CACHE_DIR, pval_preprocessing_file_name)))
        print "pval saved to file"
    # plot pvals
    pvals = [cur[1] for cur in gene_pval_pair]
    logger.info("pvals (uncorrected) below 0.05: {}".format(np.sum([True if pv < 0.05 else False for pv in pvals])))
    plot_pvalues(pvals, [i * 0.01 for i in range(101)], 0.05, "{}_01".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')]))
    plot_pvalues_log_scaled(pvals, range(100), 0.05, "{}_02".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')]))
    print "gene expression t-test done"
    #### mHGT
    print "start mHGT"
    tested_gene_list = load_gene_list(gene_list_file_name=tested_gene_list_file_name)
    total_gene_list = load_gene_list(gene_list_file_name=total_gene_list_file_name)
    lambda_group_inclusion_binary_vector = [1 if cur[0] in tested_gene_list else 0 for i, cur in
                                            enumerate(gene_pval_pair)]  # and i < true_counter
    print "gene list loading done"
    if not N:
        N = len(lambda_group_inclusion_binary_vector)
    if not B:
        B = sum(lambda_group_inclusion_binary_vector)
    HGTs = None
    non_extremer_paths_DP_table =[]
    for b in range(0,B+1):
        non_extremer_paths_DP_table.append([])
        for n in range(0, N + 1):
            non_extremer_paths_DP_table[-1].append(-1)
    # non_extremer_paths_DP_table = np.full((B+1, N+1), long(-1))
    if not hgt_preprocessing_file_name:
        hgt_preprocessing_file_name = "HGTs_out_{}_{}.npy".format(N,B)
    if os.path.isfile(os.path.join(CACHE_DIR, hgt_preprocessing_file_name)) and USE_CACHE:
        HGTs =  np.load(os.path.join(CACHE_DIR, hgt_preprocessing_file_name))
        print "hgt loaded from file"
        # left_tails = [(hypergeom.sf(i, N, B, i) + hypergeom.pmf(i, N, B, i)) for i in range(0,B+1)]
        # top_tails = [(hypergeom.sf(B, N, B, i) + hypergeom.pmf(B, N, B, i)) for i in range(0,N+1)]
    else:
        for n in range(N+1):
            # print "n1: {}".format(n)
            b_tails = np.add(hypergeom.sf(np.arange(0, B + 1), N, B, n), hypergeom.pmf(np.arange(0, B + 1), N, B, n))
            if type(HGTs)!=type(None):
                HGTs = np.c_[HGTs, b_tails]
            else:
               HGTs = b_tails
        HGTs[0][0] = 1
        np.save(os.path.join(CACHE_DIR, hgt_preprocessing_file_name), HGTs)
        print "hgt saved to file"

    print "HGTs calculations done"

    b=0
    n=0
    mHGT = 1
    for cur_step in lambda_group_inclusion_binary_vector:
       # if B == b+cur_step: break
       mHGT = min((hypergeom.sf(b, N, B, n) + hypergeom.pmf(b, N, B, n)), mHGT)
       n+=1
       if cur_step == 1:
            b=b+1
       if n > N: break

    # for i, cur in enumerate(left_tails):
    #     if cur <= mHGT:
    #         left_edge = (i,i)
    #         break
    #
    # for i, cur in enumerate(top_tails):
    #     if cur >= mHGT:
    #         top_edge = (i-1,B)
    #         break
    # slope = (left_edge[1] - top_edge[1]) / ((left_edge[0] - top_edge[0])*1.0)
    # constant = top_edge[1] - slope*top_edge[0]
    print "mHGT calculations done with {} genes included from gene family".format(b)
    for n in range(0, N+1):
        for b1 in range(0, B+1):
            # print "n: {}, b: {}".format(n, b1)
            calc_num_of_non_extremer_paths(non_extremer_paths_DP_table, HGTs, mHGT, n, b1)
         # if n !=0:
        #     for b in range(0, B + 1):
        #         non_extremer_paths_DP_table[b][n-1] = 0

    non_extremer_paths_counter = non_extremer_paths_DP_table[B][N]# calc_num_of_non_extremer_paths(non_extremer_paths_DP_table, HGTs, mHGT, N, B)

    total_possible_paths = nCk(N,B)
    print "number of valid non extremer paths: {}, total: {}, ratio: {}".format(non_extremer_paths_counter, total_possible_paths, Decimal(non_extremer_paths_counter)/Decimal(total_possible_paths))
    pval = 1- Decimal(non_extremer_paths_counter)/Decimal(total_possible_paths)
    print "mHGT pval: {}".format(pval)



    # correct pval (i.e generate qvals out of pvals)
    fdr_results = fdrcorrection0(pvals, alpha=0.05, method='indep', is_sorted=True)
    true_counter = len([cur for cur in fdr_results[0] if cur == True])
    print "true hypothesis: {}/{}".format(true_counter, np.size(fdr_results[0]))

    return pval

    # # plot qvals
    # qvals = fdr_results[1]
    # qvals.sort()
    # qval_threshold = qvals[true_counter-1]
    # plot_pvalues(qvals, [i * 0.01 for i in range(101)], qval_threshold, "{}_11".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')]))
    # plot_pvalues_log_scaled(qvals, range(100), qval_threshold, "{}_22".format(tested_gene_list_file_name[:tested_gene_list_file_name.find('.')]))

    # # summarize genes proportion
    # tested_gene_list = load_gene_list(gene_list_file_name=tested_gene_list_file_name)
    # total_gene_list = load_gene_list(gene_list_file_name=total_gene_list_file_name)
    # expected_actual_difference_list, expected_actual_ratio_difference_list, rank_in_total_list_list, \
    # z_test_proportion_test_list = summarize_genes_proportion(tested_gene_list, total_gene_list, gene_pval_pair, true_counter)


    # # calculate_gene_concentration
    # expected_tested_genes_ratio = len(tested_gene_list)*true_counter/(len(total_gene_list)*1.0)
    # plot_genes_proportion(expected_actual_difference_list, expected_actual_ratio_difference_list, z_test_proportion_test_list, rank_in_total_list_list, true_counter, expected_tested_genes_ratio, tested_gene_list_file_name)
