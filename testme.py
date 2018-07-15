
from utils.ensembl2entrez import ensembl2entrez_convertor
from matplotlib import style
import constants
style.use("ggplot")
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import os
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go

import os
import sys
import timeit
import datetime

from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.rpt.write_hierarchy import WrHierGO




def test_all():

    obo_dag = GODag(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))

    assoc = read_ncbi_gene2go(os.path.join(constants.GO_DIR, constants.GO_ASSOCIATION_FILE_NAME), no_top=True)

    """Run numerous tests for various reports."""
    dag_fin = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/mini_obo.obo")

    godag = GODag(dag_fin)
    gosubdag = GoSubDag(godag.keys(), godag)

    out = sys.stdout
    write_hier_all(gosubdag, out)