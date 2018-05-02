import re
import os
import math
from numpy import *
import numpy as np
import matplotlib as mpl
import requests
from bs4 import BeautifulSoup
from infra import *
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

    mirbase = requests.get('http://www.mirbase.org/cgi-bin/query.pl?terms={}'.format(mir_id))
    contents = mirbase.content

    mirbase_soup = BeautifulSoup(contents, 'html.parser')
    mirbase_aa = mirbase_soup.find_all('a')

    for a in mirbase_aa:
        link = a.get('href')
        if "www.targetscan.org" in link:
            taget_scan_page = requests.get(link)
            taget_scan_soup = BeautifulSoup(taget_scan_page.content, 'html.parser')
            target_scan_aa = taget_scan_soup.find_all('a')
            for a in target_scan_aa:
                if "http://www.ensembl.org/Homo_sapiens/Gene/Summary" in a.get("href"):
                    txt = a.text
                    if mirclusters.has_key(txt):
                        mirclusters[txt].add(mir_id)
                    else:
                        mirclusters[txt] = set([mir_id])
                    print txt



f = open(os.path.join(BASE_OUTPUT_DIR, "mir_cluster_out.txt"),'w')
for key, values in mirclusters.iteritems():
    line = key+"\t"
    for value in values:
        line+=(value+"\t")
    line = line[:-1] + "\n"
    f.write(line)

f.close()



