from __future__ import division
import math
import re
import sys
from itertools import product
from sklearn.externals import joblib

def load_file(file_name):
    L = []
    f = open(file_name)
    for line in f:
        if line[0] == '>':
            continue
        L.append(line.strip())
    return L

pos = load_file('data/promoterPositive.fa.txt')
neg = load_file('data/promoterNegative.fa.txt')

DUMP_DATA = False
if DUMP_DATA:
    L = ['A', 'C', 'G', 'T', '[AC]', '[AG]', '[AT]', '[CG]', '[CT]', '[GT]', '[ACGT]']
    H_pos = {}
    H_neg = {}
    for i, p in enumerate(product(L, repeat=6)):
        if i % 1000 == 0:
            print i
        patt = ''.join(p)
        key = ' '.join(p).replace(']', '').replace('[', '')
        H_pos[key] = 0
        H_neg[key] = 0
        for seq in pos:
            num_matches = len(re.findall(patt, seq))
            H_pos[key] += num_matches
        for seq in neg:
            num_matches = len(re.findall(patt, seq))
            H_neg[key] += num_matches

    joblib.dump(H_pos, 'blobs/H_pos.pkl')
    joblib.dump(H_neg, 'blobs/H_neg.pkl')
    sys.exit(0)

def get_prob_match(patt):
    p = 1
    for token in patt.split():
        p *= len(token) / 4
    return p

def get_zscore(patt, L, H):
    p_match = get_prob_match(patt)
    E_match = (len(L[0]) - 6 + 1) * len(L) * p_match
    N_match = H[patt]
    return (N_match - E_match) / math.sqrt(E_match)

H_pos = joblib.load('blobs/H_pos.pkl')
H_neg = joblib.load('blobs/H_neg.pkl')
Z = {}
for p in H_pos:
    zscore = get_zscore(p, pos, H_pos)
    Z[p] = zscore

print "          PATTERN \t POS_COUNT \t POS_ZSCORE \t NEG_COUNT \t NEG_ZSCORE"
for w in sorted(Z, key=Z.get, reverse=True)[:20]:
    print "%20s" % w, '%10.4f ' % H_pos[w], '%10.4f ' % Z[w], '%10.4f ' % H_neg[w], '%10.4f ' % get_zscore(w, neg, H_neg)

"""
        PATTERN         POS_COUNT  POS_ZSCORE  NEG_COUNT   NEG_ZSCORE
         G C G C G C     6.0000     11.9770     16.0000     16.3454
         C G C G C G     6.0000     11.9770     16.0000     16.3454
        AT A A A A A     8.0000     11.0650     11.0000      7.0864
        T G AG C G C     8.0000     11.0650      1.0000     -0.5468
       T G AG CT G C    11.0000     10.4572      2.0000     -0.7732
      T GT A A AT AT    15.0000      9.6498     35.0000     10.7378
        G C G C CG C     7.0000      9.5967     16.0000     10.9029
        AG C G C G C     7.0000      9.5967     17.0000     11.6662
        C G C G C GT     7.0000      9.5967     17.0000     11.6662
        C G C CG C G     7.0000      9.5967     16.0000     10.9029
       T T A A AT AT    10.0000      9.4190     23.0000     10.5614
       AG C G C CG C    10.0000      9.4190     17.0000      7.3229
       GT CG G C G C    10.0000      9.4190     17.0000      7.3229
     T T AT AG AT AT    21.0000      8.9749     61.0000     12.7567
      T T AT AG T AT    14.0000      8.9157     30.0000      8.8295
      A AT AT T CG T    14.0000      8.9157     22.0000      5.7763
      T GT AT A AT T    14.0000      8.9157     35.0000     10.7378
      A AT AT T GT T    14.0000      8.9157     31.0000      9.2112
      GT CG AG C G C    14.0000      8.9157     17.0000      3.8680
    T GT AT A ACGT T    20.0000      8.4558     49.0000      9.5182
"""
