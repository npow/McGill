from __future__ import division
import numpy as np
import re
from itertools import combinations
from sklearn.cross_validation import LeaveOneOut

def get_consensus(i, L):
    H = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
    for seq in L:
        H[seq[i]] += 1
    for nucleotide in H:
        if H[nucleotide] >= 0.9 * len(L):
            return nucleotide
    for pair in combinations(['A', 'C', 'T', 'G'], 2):
        if H[pair[0]] + H[pair[1]] >= 0.9 * len(L):
            return "[%s%s]" % pair
    return "[ACTG]"

def load_file(file_name):
    L = []
    f = open(file_name)
    for line in f:
        L.append(line.strip())
    return L

def get_patt(L):
    patt = ''
    for i in range(6):
        patt += get_consensus(i, L)
    return patt

def main(dataset):
    print "dataset%d" % dataset
    pos = np.array(load_file('data/dataset%d/positive.txt' % dataset))
    neg = np.array(load_file('data/dataset%d/negative.txt' % dataset))

    TP = 0
    TN = 0
    loo = LeaveOneOut(len(pos))
    for train, test in loo:
        patt = get_patt(pos[train])
        if re.match(patt, pos[test][0]):
            TP += 1
    patt = get_patt(pos)
    for seq in neg:
        if re.match(patt, seq) is None:
            TN += 1
    print "sensitivity: %f" % (TP/len(pos))
    print "specificity: %f" % (TN/len(neg))

if __name__ == '__main__':
    main(1)
    main(2)
