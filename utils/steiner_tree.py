"""Print a GO term's lower-level hierarchy."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx
from goatools.godag.consts import Consts
from goatools.gosubdag.go_paths import GoPaths
import constants
"""driver imports"""
import os
import timeit
import numpy as np
import pcst_fast
import datetime
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.rpt.write_hierarchy import WrHierGO
""""gene to go mapper imports"""
import time
from goatools import obo_parser
from  goatools.associations import read_ncbi_gene2go
from goatools.base import download_ncbi_associations
import constants
from utils.ensembl2entrez import get_entrez2ensembl_dictionary
import wget, os
from utils.ensg_dictionary import get_ensg_dict
from utils.ensp_dictionary import get_ensp_dict
from utils.string_ppi_dictionary import get_string_ppi_dict
from neo4j.v1 import GraphDatabase
import infra

class WrHierGO(object):

    def __init__(self, gosubdag, **kws):
        self.gosubdag = gosubdag  # GoSubDag arg, children=True, must be used
        self.usrdct = {k:v for k, v in kws.items() if k in kws}
        self.usrset = set([k for k, v in kws.items() if k in kws and v])
        # ' {NS} {dcnt:6,} L{level:02} D{depth:02} {D1:5} {GO_name}'


    def ext_hier_down(self, goid, prt=sys.stdout):

        obj = _WrHierPrt(self, prt)
        obj.ext_hier_rec(goid)
        prt.write("VERTICES:\n")
        for k, v in obj.vertices.iteritems():
            prt.write(str("{} - {}\n".format(k,v)))

        prt.write("EDGES:\n")
        for k, v in obj.edges.iteritems():
            prt.write(str("{} - {}\n".format(k,v)))

        return obj

class _WrHierPrt(object):

    def __init__(self, obj, prt=sys.stdout):
        self.gosubdag = obj.gosubdag
        self.max_indent = obj.usrdct.get('max_indent')
        self.include_only = obj.usrdct['include_only'] if 'include_only' in obj.usrdct else None
        self.go_marks = obj.usrdct['go_marks'] if 'go_marks' in obj.usrdct else set()
        self.concise_prt = 'concise' in obj.usrset
        self.indent = 'no_indent' not in obj.usrset
        self.go2geneids = obj.usrdct.get('go2geneids')
        # vars
        self.prt = prt
        self.edges = {}
        self.vertices = {}
        self.gos_printed = set()
        self.prtfmt = self._init_prtfmt()
        self.dash_len = obj.usrdct.get('dash_len', 6) + 12

    def ext_hier_rec(self, goid, depth=1):
        """Write hierarchy for a GO Term record and all GO IDs down to the leaf level."""
        ntgo = self.gosubdag.go2nt[goid]
        ntobj = self.gosubdag.go2obj[goid]
        # Shortens hierarchy report by only printing the hierarchy
        # for the sub-set of user-specified GO terms which are connected.
        if self.include_only and goid not in self.include_only:
            return
        nrp = self.concise_prt and goid in self.gos_printed
        if self.vertices.has_key(goid):
            self.vertices[goid]["weight"] += 1
            self.vertices[goid]["depth"].append(depth)
        else:
            self.vertices[goid] = {"name": ntgo.GO_name,
                                   "weight": 0,
                                   "NS": ntgo.NS,
                                   "depth": [depth],
                                   "L" : ntgo.level,
                                   "D" : ntgo.depth,
                                   "isleaf" : len(ntobj.children) == 0,
                                   "obj": ntobj}

        self.gos_printed.add(goid)
        # Do not extract hierarchy below this turn if it has already been printed
        # if nrp:
        #     return
        depth += 1
        if self.max_indent is not None and depth > self.max_indent:
            return
        for child in ntobj.children:
            if self.go2geneids.has_key(child.id) or True:
                if self.edges.has_key("{}={}".format(goid, child.id)):
                    self.edges["{}={}".format(goid, child.id)]["weight"] += 1
                else:
                    self.edges["{}={}".format(goid, child.id)] = {"weight" : 0}
                self.ext_hier_rec(child.id, depth)


def extract_hier_all(gosubdag, out, root_term, go2geneids):
    """write_hier.py: Prints the entire mini GO hierarchy, with counts of children."""
    out.write('\nTEST EXTRACTION: Print all hierarchies:\n')
    objwr = WrHierGO(gosubdag,  go2geneids=go2geneids)
    obj = objwr.ext_hier_down(root_term, out)
    return (obj.vertices, obj.edges)

def fetch_go_hierarcy():
    obo_file_location = os.path.join(constants.GO_DIR, constants.GO_FILE_NAME)
    if not os.path.exists(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME)):
        wget.download(constants.GO_OBO_URL, os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))

    print "Downloading gene-GO associations"
    association_file_location = os.path.join(constants.GO_DIR, constants.GO_ASSOCIATION_FILE_NAME)
    if not os.path.exists(association_file_location):
        wget.download(constants.GO_ASSOCIATION_GENE2GEO_URL,
                      os.path.join(constants.GO_DIR, constants.GO_ASSOCIATION_FILE_NAME))

    print "Loading gene-GO associations"
    # gene2go = download_ncbi_associations(obo_file_location) - why does this line needed?
    go2geneids = read_ncbi_gene2go(association_file_location, taxids=[9606], go2geneids=True)
    geneids2go = read_ncbi_gene2go(association_file_location, taxids=[9606])

    return (go2geneids, geneids2go)

def fetch_string_ppi_edges():
    go_edges = {}
    if constants.USE_CACHE:
        if os.path.isfile(os.path.join(constants.DICTIONARIES_DIR, "GO_edges_ppi_filtered_cc_leafs.txt")):
            GO_edges_ppi_grid = infra.load_phenotype_data("GO_edges_ppi_filtered_cc_leafs.txt",phenotype_list_path=constants.DICTIONARIES_DIR)
            for cur in GO_edges_ppi_grid:
                go_edges[cur[0]] = int(cur[1])
            return go_edges

    print "fetching ensg"
    ensg_dict = get_ensg_dict()
    print "fetching ensp"
    ensp_dict = get_ensp_dict()
    print "fetching string ppi"
    string_ppi_dict = get_string_ppi_dict()
    go_edges = {}
    count = 0
    for cur_edge, cur_score in string_ppi_dict.iteritems():
        count +=1
        print count
        vertices = cur_edge.split("=")
        if not ensp_dict.has_key(vertices[0]) or not ensp_dict.has_key(vertices[1]): continue

        go_src = ensp_dict[vertices[0]]["GO Terms"]
        go_dst = ensp_dict[vertices[1]]["GO Terms"]

        for cur_src in go_src:
            for cur_dst in go_dst:
                edge = "{}={}".format(cur_src,cur_dst)
                edge_alt = "{}={}".format(cur_dst, cur_src)
                if go_edges.has_key(edge):
                    go_edges[edge] += int(cur_score)
                elif go_edges.has_key(edge_alt):
                    go_edges[edge_alt] += int(cur_score)
                else:
                    go_edges[edge] = int(cur_score)
    with file(os.path.join(constants.OUTPUT_GLOBAL_DIR, "GO_edges_ppi_filtered_cc_leafs.txt"), "w+") as f:
        for k,v in go_edges.iteritems():
            f.write("{}\t{}\n".format(k,v))

    return go_edges



#################################################################
# Driver
#################################################################
def build_hierarcy():
    print "fetching ppi"
    go_edges = fetch_string_ppi_edges()

    go2geneids, geneids2go = fetch_go_hierarcy()

    """Run numerous tests for various reports."""
    dag_fin = os.path.join(constants.GO_DIR, constants.GO_FILE_NAME)
    tic = timeit.default_timer()
    godag = GODag(dag_fin, optional_attrs=['relationship'])
    gosubdag = GoSubDag(godag.keys(), godag)
    toc = timeit.default_timer()
    out = file(os.path.join(constants.BASE_PROFILE, "output", "go_hierarcy.txt"), "w+")  # sys.stdout
    dict_result = {}
    for cur_term in ['GO:0005575']:
        vertices, edges = extract_hier_all(gosubdag, out, cur_term, go2geneids)
        dict_result[cur_term] = {"vertices": vertices, "edges": edges}



    def add_edge(tx, src, dst, score):
        tx.run(("MERGE (n1: GO{{term:\"{TERM1}\"}})"+ \
                "MERGE (n2: GO{{term:\"{TERM2}\"}})"+ \
                "MERGE (n1)-[r:SCR {{ score: {SCORE} }}]->(n2)").format(TERM1=src, TERM2=dst, SCORE=score))

    def add_node(tx, nd):
        tx.run(("CREATE (n1: GO{{term:\"{TERM1}\"}})".format(TERM1=nd)))


        count=0
        vertices_grid = []
        vertices_prizes = []
        for k, v in dict_result['GO:0005575']['vertices'].iteritems():
            if dict_result['GO:0005575']['vertices'].has_key(k) \
                            and dict_result['GO:0005575']['vertices'][k]['isleaf']:
                vertices_grid.append(k)
                vertices_prizes.append(1)
                count+=1
        print "total vartices: {}".format(count)


        count=0
        edges_grid = []
        edges_costs = []
        for cur_edges, score in go_edges.iteritems():

            vertices = cur_edges.split("=")
            if dict_result['GO:0005575']['vertices'].has_key(vertices[0]) and dict_result['GO:0005575'][
                'vertices'].has_key(vertices[1]) and score > 100000 \
                    and dict_result['GO:0005575']['vertices'][vertices[0]]['isleaf'] and \
                    dict_result['GO:0005575']['vertices'][vertices[1]]['isleaf']:
                edges_grid.append([vertices_grid.index(vertices[0]), vertices_grid.index(vertices[0])])
                edges_costs.append(1/float(score))
                count+=1
        print "total edges: {}".format(count)

        edges_grid=np.array(edges_grid).astype(np.float64)
        vertices_prizes=np.array(vertices_prizes).astype(np.floa64)
        edges_costs = np.array(edges_costs).astype(np.float64)
        root = -1
        num_clusters=1
        pruning = None
        verbosity_level = None
        vertices, edges = pcst_fast(edges_grid, vertices_prizes, edges_costs, root, num_clusters, pruning, verbosity_level)



if __name__ == '__main__':
    print build_hierarcy()

