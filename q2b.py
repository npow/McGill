from __future__ import division
import numpy as np
import re
from itertools import combinations
from sklearn.cross_validation import LeaveOneOut

def get_pwm(L):
    H = {'A': [0]*6, 'C': [0]*6, 'T': [0]*6, 'G': [0]*6}
    for i in xrange(6):
        for seq in L:
            H[seq[i]][i] += 1
    for h in H:
        for i in xrange(6):
            H[h][i] /= len(L)
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

def main(dataset):
    print "dataset%d" % dataset
    pos = np.array(load_file('data/dataset%d/positive.txt' % dataset))
    neg = np.array(load_file('data/dataset%d/negative.txt' % dataset))

    loo = LeaveOneOut(len(pos))
    pos_proba = []
    for train, test in loo:
        H = get_pwm(pos[train])
        p = f(H, pos[test][0])
        pos_proba.append(p)
    H = get_pwm(pos)
    neg_proba = []
    for seq in neg:
        p = f(H, seq)
        neg_proba.append(p)

    min_errors = float('inf')
    best = (0, 0)
    best_thresh = None
    for t in xrange(0, 10001):
        threshold = t / 10000.0
        TP = len(filter(lambda x: x > threshold, pos_proba))
        TN = len(filter(lambda x: x <= threshold, neg_proba))
        FP = len(pos) - TP
        FN = len(neg) - TN
        if FP + FN < min_errors:
            min_errors = FP + FN
            best = (TP, TN)
            best_thresh = threshold

    print "threshold: %f" % best_thresh
    print "sensitivity: %f" % (best[0]/len(pos))
    print "specificity: %f" % (best[1]/len(neg))

if __name__ == '__main__':
    main(1)
    main(2)
