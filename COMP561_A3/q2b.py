from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import re
from itertools import combinations
from sklearn.cross_validation import LeaveOneOut
from sklearn import metrics

def get_pwm(L, ALPHA):
    H = {'A': [ALPHA]*6, 'C': [ALPHA]*6, 'T': [ALPHA]*6, 'G': [ALPHA]*6}
    for i in xrange(6):
        for seq in L:
            H[seq[i]][i] += 1
    for h in H:
        for i in xrange(6):
            H[h][i] /= len(L) + 4*ALPHA
    return H

def f(H, seq):
    p = 1
    for i in range(6):
        p *= H[seq[i]][i]
    return p

def load_file(file_name):
    L = []
    f = open(file_name)
    for line in f:
        L.append(line.strip())
    return L

def main(dataset, alpha):
    pos = np.array(load_file('data/dataset%d/positive.txt' % dataset))
    neg = np.array(load_file('data/dataset%d/negative.txt' % dataset))

    loo = LeaveOneOut(len(pos))
    pos_proba = []
    for train, test in loo:
        H = get_pwm(pos[train], alpha)
        p = f(H, pos[test][0])
        pos_proba.append(p)
    H = get_pwm(pos, alpha)
    neg_proba = []
    for seq in neg:
        p = f(H, seq)
        neg_proba.append(p)

    min_errors = float('inf')
    best = (0, 0)
    best_thresh = None
    thresholds = []
    sensitivities = []
    specificities = []
    for threshold in sorted([0, 1] + pos_proba + neg_proba):
        TP = len(filter(lambda x: x >= threshold, pos_proba))
        TN = len(filter(lambda x: x < threshold, neg_proba))
        FP = len(pos) - TP
        FN = len(neg) - TN
        thresholds.append(threshold)
        sensitivities.append(TP/len(pos))
        specificities.append(TN/len(neg))
        if FP + FN < min_errors:
            min_errors = FP + FN
            best = (TP, TN)
            best_thresh = threshold

    print "*" * 80
    print "dataset%d" % dataset
    print "alpha: %d" % alpha
    print "threshold: %f" % best_thresh
    print "sensitivity: %f" % (best[0]/len(pos))
    print "specificity: %f" % (best[1]/len(neg))
    auc = metrics.auc(sensitivities, specificities, reorder=True)
    print "auc: %f" % auc
    plt.plot(sensitivities, specificities, linewidth=4)
    #plt.xscale('log')
    plt.xlabel('Sensitivity')
    plt.ylabel('Specificity')
    plt.title('Sensitivity-Specificity, dataset%d. AUC: %f' % (dataset, auc))
    plt.legend()
    plt.savefig('plots/ss%s_dataset%d' % ('' if alpha == 0 else '_smooth', dataset))
    plt.clf()

if __name__ == '__main__':
    main(1, alpha=0)
    main(2, alpha=0)
    main(1, alpha=1)
    main(2, alpha=1)
