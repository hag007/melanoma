import os
import sys
import timeit
import datetime

from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.rpt.write_hierarchy import WrHierGO

#################################################################
# Sub-routines to tests
#################################################################
def write_hier_all(gosubdag, out):
    """write_hier.py: Prints the entire mini GO hierarchy, with counts of children."""
    out.write('\nTEST ALL: Print all hierarchies:\n')
    objwr = WrHierGO(gosubdag)
    gos_printed = objwr.prt_hier_down("GO:0000001", out)
    assert gos_printed == set(objwr.gosubdag.go2nt)


def write_hier_norep(gosubdag, out):
    """Shortens hierarchy report by only printing branches once.
         Prints the 'entire hierarchy' of GO:0000005 the 1st time seen:
           ---     1 GO:0000005    L-02    D-02
           ----     0 GO:0000010   L-03    D-04
         Prints just GO:0000005 (ommit child GO:10) the 2nd time seen:
           ===     1 GO:0000005    L-02    D-02
         '=' is used in hierarchy mark to indicate that the pathes
             below the marked term have already been printed.
    """
    out.write('\nTEST ALL: Print branches just once:\n')
    objwr = WrHierGO(gosubdag, concise=True)
    gos_printed = objwr.prt_hier_down("GO:0000001", out)
    assert gos_printed == set(objwr.gosubdag.go2nt)


def write_hier_lim(gosubdag, out):
    """Limits hierarchy list to GO Terms specified by user."""
    go_omit = ['GO:0000005', 'GO:0000010']
    go_ids = [go_id for go_id in gosubdag.go2obj if go_id not in go_omit]
    out.write('\nTEST OMIT: 05 and 10:\n')
    objwr = WrHierGO(gosubdag, include_only=go_ids)
    gos_printed = objwr.prt_hier_down("GO:0000001", out)
    assert not gos_printed.intersection(go_omit), "SHOULD NOT PRINT {GOs}".format(GOs=go_omit)


def write_hier_mrk(gosubdag, out):
    """Print all paths, but mark GO Terms of interest. """
    mark_lst = ['GO:0000001', 'GO:0000003', 'GO:0000006', 'GO:0000008', 'GO:0000009']
    out.write('\nTEST MARK: 01->03->06->08->09:\n')
    objwr = WrHierGO(gosubdag, go_marks=mark_lst)
    objwr.prt_hier_down("GO:0000001", out)
      #go_marks=[oGO.id for oGO in oGOs_in_cluster])

#################################################################
# Tests
#################################################################
def test_all():
    """Run numerous tests for various reports."""
    dag_fin = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/mini_obo.obo")
    tic = timeit.default_timer()
    godag = GODag(dag_fin)
    gosubdag = GoSubDag(godag.keys(), godag)
    toc = timeit.default_timer()
    out = sys.stdout
    write_hier_all(gosubdag, out)
    write_hier_norep(gosubdag, out)
    write_hier_lim(gosubdag, out)
    write_hier_mrk(gosubdag, out)
    msg = "Elapsed HMS: {}\n\n".format(str(datetime.timedelta(seconds=(toc-tic))))
    sys.stdout.write(msg)

#################################################################
# main
#################################################################
if __name__ == '__main__':
    test_all()
