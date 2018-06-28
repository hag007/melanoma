import Bio.UniProt.GOA as GOA
# from orangecontrib.bio.go import Ontology
import wget
from utils.ensembl2gene_symbol import e2g_convertor
import time
import requests
import scipy.special
import matplotlib.pyplot as plt
from matplotlib import style
import matplotlib.ticker as ticker
style.use("ggplot")
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from infra import *
import os
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import proj3d

def plot_pca(dataset, labels_assignment, meta_groups):
    actual_labels = list(range(1,len(meta_groups[0])+1))

    labels = [(cur["_name"], cur["_label"]) for i, cur in enumerate(meta_groups[0])]
    labels = [("unknown", 0)]+labels
    X=[]
    y=[]
    for cur in actual_labels:
        X.append(np.average(dataset[np.where(cur == labels_assignment[0]), :][0], axis=0))
        y.append(cur)
    X = np.array(X)
    y = np.array(y)
    pca = PCA(n_components=3)
    pca.fit_transform(X[:len(X)-1])

    fig = plt.figure(1, figsize=(20, 20))
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='jet')
    for i, x in enumerate(X):
        name = labels[i+1]
        x2, y2, _ = proj3d.proj_transform(x[0], x[1], x[2], ax.get_proj())
        ax.annotate(name,
                    xy=(x2, y2), xytext=(-20, 20), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))

    plt.savefig(os.path.join(constants.BASE_PROFILE, "output", "PCA_{}_{}.png").format(constants.CANCER_TYPE,time.time()))

def plot_pca_by_samples(dataset, labels_assignment, meta_groups):
    # actual_labels = list(range(1,len(meta_groups[0])+1))

    # labels = [(cur["_name"], cur["_label"]) for i, cur in enumerate(meta_groups[0])]
    # labels = [("unknown", 0)]+labels
    # X=[]
    # y=[]
    # for cur in actual_labels:
    #     X.append(np.average(dataset[np.where(cur == labels_assignment[0]), :][0], axis=0))
    #     y.append(cur)
    # X = np.array(X)
    # y = np.array(y)
    pca = PCA(n_components=3)
    pca.fit_transform(dataset)

    fig = plt.figure(1, figsize=(15, 15))
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(dataset[:, 0], dataset[:, 1], dataset[:, 2])
    # for i, x in enumerate(X):
    #     name = labels[i+1]
    #     x2, y2, _ = proj3d.proj_transform(x[0], x[1], x[2], ax.get_proj())
    #     ax.annotate(name,
    #                 xy=(x2, y2), xytext=(-20, 20), textcoords='offset points',
    #                 bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
    #                 arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=3, headlength=2))

    plt.savefig(os.path.join(constants.BASE_PROFILE, "output", "PCA_{}_{}.png").format(constants.CANCER_TYPE,time.time()))