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

pvals = [0.9813, 0.4420, 0.2890, 0.9791, 0.9991, 0.7856, 0.0314, 0.6535, 0.1499, 0.2626, 0.0437]


pvals.sort()
fdr_results = fdrcorrection0(pvals, alpha=0.05, method='indep', is_sorted=True)
true_counter = len([cur for cur in fdr_results[0] if cur == True])
print "true hypothesis: {}".format(true_counter)
print "total hypothesis: {}".format(np.size(fdr_results[0]))