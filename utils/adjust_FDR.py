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

# pvals = [0.9813, 0.4420, 0.2890, 0.9791, 0.9991, 0.7856, 0.0314, 0.6535, 0.1499, 0.2626, 0.0437]
#counts
# pvals = [2.72E-01, 9.73E-01, 3.27E-06, 3.79E-02, 2.18E-01, 9.09E-01]
# fpkm
pvals = [0.122144119919258, 0.195818023252933, 1.07642059824493E-07, 0.0124414358766162, 0.0144899123857276, 0.611362249862538]
# fpkm-uq
# pvals = [0.04364191582855, 0.0746803276988296, 3.23745624934477E-09, 0.0337253006483806, 0.00202830018481174, 0.0278338429119658]
pvals.sort()
fdr_results = fdrcorrection0(pvals, alpha=0.05, method='indep', is_sorted=True)
true_counter = len([cur for cur in fdr_results[0] if cur == True])
print fdr_results[1]
print "true hypothesis: {}".format(true_counter)
print "total hypothesis: {}".format(np.size(fdr_results[0]))