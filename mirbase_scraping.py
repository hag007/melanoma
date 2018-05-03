import re
import os
import math
from numpy import *
import numpy as np
import matplotlib as mpl
import requests
from bs4 import BeautifulSoup
from infra import *

def load_gene_dictionary(gene_list_file_name, gene_list_path=None, source="GDC-TCGA",dataset="melanoma"): #  ="TCGA-SKCM.htseq_counts.tsv"
    if gene_list_path == None:
        gene_list_path = os.path.join(BASE_PROFILE,source,dataset,gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines

def send_request(link):
    counter = 0
    req_success = False
    while not req_success:
        try:
            page = requests.get(link)
            soup = BeautifulSoup(page.content, 'html.parser')
            req_success = True
        except Exception:
            counter+=1
            print "Failed to establish a new connection to {}. retry {} time.".format(link, counter)
            pass
    return soup


lines_dict = load_gene_dictionary("ensembl2gene_symbol.txt")


gene_symbols2ensembl = {}
included_genes = []
for cur in lines_dict:
    splited_line = cur.split()
    if splited_line[0].find('.') > 0:
        limit = splited_line[0].find('.')
    else:
        limit = len(splited_line[0])
    gene_symbols2ensembl[splited_line[1]] = splited_line[0][:limit]

mirclusters = {}

mir_ids = load_gene_list("mir_total.txt")
for mir_id in mir_ids:

    mirbase_soup = send_request('http://www.mirbase.org/cgi-bin/query.pl?terms={}'.format(mir_id))
    mirbase_aa = mirbase_soup.find_all('a')

    for a in mirbase_aa:
        link = a.get('href')
        if "www.targetscan.org" in link:
            taget_scan_soup = send_request(link)
            target_scan_aa = taget_scan_soup.find_all('a')
            for a in target_scan_aa:
                if "http://www.ensembl.org/Homo_sapiens/Gene/Summary" in a.get("href"):
                    txt = a.text
                    if mirclusters.has_key(txt):
                        mirclusters[txt].add(mir_id)
                    else:
                        mirclusters[txt] = set([mir_id])
                    # print txt
    print "done analyzing {}".format(mir_id)


f = open(os.path.join(BASE_OUTPUT_DIR, "mir_cluster_out.txt"),'w')
for key, values in mirclusters.iteritems():
    if gene_symbols2ensembl.has_key(key):
        line = gene_symbols2ensembl[key]+"\t"
        for value in values:
            line+=(value+"\t")
        line = line[:-1] + "\n"
        f.write(line)

f.close()



