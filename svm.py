import sys
import mord
import scipy.special
from matplotlib import style
style.use("ggplot")
from sklearn import svm
from sklearn.model_selection import GridSearchCV, cross_val_score
import scipy
from scipy.stats import hypergeom
import logging
sh = logging.StreamHandler()
logger = logging.getLogger("log")
logger.addHandler(sh)
from constants import *
from infra import *
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve


DISTANCE = "distance"
LOGISTIC_REGRESSION = "logistic_regression"


############################################(3) SVM to see whether we can predict tumor type by tested genes' profile expression########################################

def load_svm_data(tested_gene_file_name, expression_profile_file_name, phenotype_file_name, label=None, label_values=None, gene_filter_file_name=None, groups=None):
    data = []
    labels = []
    labeled_groups = load_expression_profile_by_labelling(tested_gene_file_name, expression_profile_file_name, phenotype_file_name, label, label_values, gene_filter_file_name, groups=groups)
    print "groups sizes: {}".format([len(cur_group) for cur_group in labeled_groups])
    labeled_groups = sorted(labeled_groups, key = lambda x: len(x))

    for i, cur_group in enumerate(labeled_groups):
        for cur in cur_group[1:]:
            data.append(cur[1:].astype(np.float))
            labels.append(i)


    return data, labels, labeled_groups, labeled_groups[0][0][1:]

def randonize_patients(data, labels):

    while True:
        total  = zip(data, labels)
        random.shuffle(total)
        data, label = zip(*total)
        if sum(label[(3 * len(label)) / 4:]) % len(label[(3 * len(label)) / 4:]) != 0 and \
            sum(label[:(3 * len(label)) / 4]) % len(label[:(3 * len(label)) / 4]) != 0:
            break
    return (data, label)
def randonize_patients_old(data, labels):
    total = np.c_[data, labels]
    random.shuffle(total)
    # results = zip(*total)
    return total[:,:-1], total[:,-1]


def divide_train_and_test_groups(data, labels):
    data_train = data[:(3 * len(data)) / 4]
    data_test = data[(3 * len(data)) / 4:]
    labels_train = labels[:(3 * len(data)) / 4]
    labels_test = labels[(3 * len(data)) / 4:]
    return data_train, data_test, labels_train, labels_test

def svm_rbf_default(tuned_parameters):
    return GridSearchCV(svm.SVC(probability=True), param_grid=tuned_parameters, return_train_score=True)

def svm_multiclass(tuned_parameters):
    return mord.LogisticAT(alpha=1.)

def apply_svm(clf_method, data_train, labels_train, data_test, labels_test, rank_method):

    data_train = [[cur2 for cur2 in cur1] for cur1 in data_train]
    labels_train = [cur for i, cur in enumerate(labels_train)]
    # alternated: labels_train = [i%2 for i, cur in enumerate(labels_train)]
    # inverted: labels_train = [abs(cur - 1) for i, cur in enumerate(labels_train)]
    # randomized: labels_train = [floor(random.random()/0.5) for i, cur in enumerate(labels_train)]

    data_test = [[cur2 for cur2 in cur1] for cur1 in data_test]
    labels_test = [cur for i, cur in enumerate(labels_test)]
    # alternated: labels_test = [i%2 for i, cur in enumerate(labels_test)]
    # inverted: labels_test = [abs(cur-1) for i, cur in enumerate(labels_test)]
    # randomized: labels_test = [floor(random.random()/0.5) for i, cur in enumerate(labels_test)]

    # shuffled: random.shuffle(labels_train)
    # shuffled: random.shuffle(labels_test)

    data_train=np.array(data_train)
    labels_train=np.array(labels_train)
    data_test = np.array(data_test)
    labels_test = np.array(labels_test)
    clf_method.fit(data_train, labels_train)
    predicted_results = clf_method.predict(data_test)
    if rank_method == DISTANCE:
        probabilities = clf_method.decision_function(data_test)
    else:
        probabilities = clf_method.predict_proba(data_test)
        probabilities = probabilities[:,1]
    #zipped = zip(probabilities, data_train, labels_train)
    #zipped_sorted = sorted(zipped, key=lambda x: x[0])
    #patients_rank_sorted = [x[0] for x in zipped_sorted]
    #patients_expression_sorted = [x[1] for x in zipped_sorted]
    #patients_labels_sorted = [x[2] for x in zipped_sorted]
    precision, recall, _ = precision_recall_curve(labels_test, probabilities)
    average_precision = average_precision_score(labels_test, probabilities)
    fpr, tpr, _ = roc_curve(labels_test, probabilities)
    auc = roc_auc_score(labels_test, probabilities)


    # prediction = clf_method.predict(data_test)
    # predicted_positives = sum(prediction)
    # predicted_negatives = len(prediction) - predicted_positives
    # labeled_positives = sum(labels_test)
    # labeled_negatives = len(labels_test) - labeled_positives
    #
    # tp = []
    # fp = []
    # fn = []
    # tn = []
    # precision_1 = []
    # recall_1 = []
    # for ind, cur in enumerate(prediction):
    #     labeled_positives = sum(labels_test[:ind])
    #     labeled_negatives = len(labels_test[:ind]) - labeled_positives
    #     tp.append(sum([1 for i, cur in enumerate(prediction[:ind]) if cur == 1 and labels_test[i] == 1]))
    #     fp.append(sum([1 for i, cur in enumerate(prediction[:ind]) if cur == 1 and labels_test[i] == 0]))
    #     fn.append(labeled_positives - tp[-1])
    #     tn.append(labeled_negatives - fp[-1])
    #     precision_1.append(float(tp[-1]) / (tp[-1] + fn[-1]))
    #     recall_1.append(float(tp[-1]) / (tp[-1] + fp[-1]))
    #
    # recall_1, precision_1 = zip(*sorted(zip(recall, precision), key=lambda x: x[0]))
    #
    # accuracy = 1 - float(sum(abs(prediction - labels_test)))/len(labels_test)

    print "PR: {}, ROC: {}".format(average_precision,auc) ## , fp: {}, fn: {}, total: {}  ## , fp[-1], fn[-1], len(prediction))

    # plt.plot(recall_1, precision_1, label="PRAUC")
    # plt.step(recall, precision, color='b', alpha=0.2,
    #          where='post')
    # plt.fill_between(recall, precision, step='post', alpha=0.3,
    #                  color='b')
    #
    # plt.xlabel('Recall')
    # plt.ylabel('Precision')
    # plt.ylim([0.0, 1.05])
    # plt.xlim([0.0, 1.0])
    # plt.title('2-class Precision-Recall curve: AP={0:0.2f}'.format(average_precision))
    #
    # plt.show()
    # plt.savefig(os.path.join(BASE_OUTPUT_DIR,"test.png"))
    #################
    # plt.plot(fpr, tpr, label="PRAUC")
    # plt.step(fpr, tpr, color='b', alpha=0.2,
    #          where='post')
    # plt.fill_between(fpr, tpr, step='post', alpha=0.3,
    #                  color='b')
    #
    # plt.xlabel('fpr')
    # plt.ylabel('tpr')
    # plt.ylim([0.0, 1.05])
    # plt.xlim([0.0, 1.0])
    # plt.title('2-class Precision-Recall curve: AP={0:0.2f}'.format(auc))
    #
    # plt.show()
    # plt.savefig(os.path.join(BASE_OUTPUT_DIR,"test.png"))


    #################
    ##
    # ta = (len(labels_test) * 1.0 - sum([abs(p - r) for p, r in zip(predicted_results, labels_test)])) / (len(labels_test) * 1.0)
    # return clf_method.best_score_, ta
    return average_precision, auc

def print_svm_results(train_scores, train_alt_scores, test_scores, test_alt_scores, rounds):
    print "train_accuracy_avg: {}".format(sum(train_scores)/ rounds)
    print "train_alt_accuracy_avg: {}".format(sum(train_alt_scores)/ rounds)
    print "train_accuracy_diff_avg: {}".format((sum(train_scores)-sum(train_alt_scores))/rounds)
    print "test_accuracy_avg: {}".format(sum(test_scores)/ rounds)
    print "test_alt_accuracy_avg: {}".format(sum(test_alt_scores)/ rounds)
    print "test_accuracy_diff_avg: {}".format((sum(test_scores)-sum(test_alt_scores))/rounds)
    print "train p val: {}".format(scipy.stats.ttest_ind(train_scores, train_alt_scores)[1])
    print "test_p_val: {}".format(scipy.stats.ttest_ind(test_scores, test_alt_scores)[1])

    # print "{}".format(train_score_sum / rounds)
    # print "{}".format(train_alt_score_sum / rounds)
    # print "{}".format(scipy.stats.ttest_ind(train_scores, train_alt_scores)[1])
    # print "{}".format(test_score_sum / rounds)
    # print "{}".format(test_alt_score_sum / rounds)
    # print "{}".format(scipy.stats.ttest_ind(test_scores, train_alt_scores)[1])

def print_svm_results(test_scores, rounds):
    avgs = []
    print "start avg and var:"
    for i, cur_i in enumerate(test_scores):
        avg = sum(test_scores[i]) / rounds
        print "{}\t {}".format(avg, np.var(test_scores[i]))
        avgs.append(avg)
    print "done PR-AUC: avg and var"

    results_summarized = []
    for i, cur_i in enumerate(avgs):
        results_summarized.append([])

    for i, cur_i in enumerate(avgs):
        for j, cur_j in enumerate(avgs):
            if j>i:
                pval = scipy.stats.ttest_ind(test_scores[i], test_scores[j])[1]
                #print "test_p_val for {}, {}: {}".format(i,j,pval)
                results_summarized[i].append(pval)
            elif j==i:
                results_summarized[i].append(1)
            else:
                results_summarized[i].append(-1)

    for i, cur_i in enumerate(avgs):
        for j, cur_j in enumerate(avgs):
            if results_summarized[i][j] == -1:
                results_summarized[i][j] = results_summarized[j][i]

    # for i, cur_i in enumerate(avgs):
    #     print ""
    #     for j, cur_j in enumerate(avgs):
    #         print "{}\t".format(results_summarized[i][j]),

    # print "{}".format(train_score_sum / rounds)
    # print "{}".format(train_alt_score_sum / rounds)
    # print "{}".format(scipy.stats.ttest_ind(train_scores, train_alt_scores)[1])
    # print "{}".format(test_score_sum / rounds)
    # print "{}".format(test_alt_score_sum / rounds)
    # print "{}".format(scipy.stats.ttest_ind(test_scores, train_alt_scores)[1])


# (3) main
def prediction_by_gene_expression(tested_gene_file_names, expression_profile_file_name, phenotype_file_name, label=None, label_values=None, rank_method=LOGISTIC_REGRESSION, gene_filter_file_name=None, rounds=2, groups=None, classification_method="svm_rbf_default", tuning_parameters={'C': [10], 'kernel': ['rbf']}):
    thismodule = sys.modules[__name__]
    clf = getattr(thismodule, classification_method)(tuning_parameters)
    genelist_datasets = []
    for tested_gene_file_name in tested_gene_file_names:
        data, labels , _1, _2 = load_svm_data(tested_gene_file_name, expression_profile_file_name, phenotype_file_name, label, label_values,
                                     gene_filter_file_name, groups)
        genelist_datasets.append(data)

    train_scores = []
    test_pr_score = []
    test_roc_score = []
    for j in range(len(genelist_datasets)):
        train_scores.append([])
        test_pr_score.append([])
        test_roc_score.append([])
    for i in range(rounds):
        genelist_datasets = np.rot90(genelist_datasets, k=1, axes=(1, 0))
        genelist_datasets, labels = randonize_patients(genelist_datasets, labels)
        genelist_datasets = np.rot90(genelist_datasets, k=-1, axes=(1, 0))
        for j, cur_dataset in enumerate(genelist_datasets):
            data_train, data_test, labels_train, labels_test = divide_train_and_test_groups(cur_dataset, labels)

            test_pr, test_roc = apply_svm(clf, data_train, labels_train, data_test, labels_test, rank_method)
            test_pr_score[j].append(test_pr)
            test_roc_score[j].append(test_roc)
    print "#######################################"
    print "PRAUC"
    print_svm_results(test_pr_score, float(rounds))
    print "ROC"
    print_svm_results(test_roc_score, float(rounds))
