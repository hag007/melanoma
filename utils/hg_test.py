from scipy.stats import hypergeom

def calc_HG_test(total_gene_list_N, tests_gene_list_B, total_gene_list_n):
    b = len(set(total_gene_list_n).intersection(tests_gene_list_B))
    B = len(tests_gene_list_B)
    N = len(total_gene_list_N)
    n = len(total_gene_list_n)
    print "run HG test with {},{},{},{}".format(b, N, B, n)
    return hypergeom.sf(b - 1, N, B, n)