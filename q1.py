from __future__ import division
import math
import re
import sys
from itertools import permutations
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
    for i, p in enumerate(permutations(L, 6)):
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

for w in sorted(Z, key=Z.get, reverse=True)[:20]:
  print "%20s" % w, '\t\t', H_pos[w], '\t', "%0.4f" % Z[w], '\t', H_neg[w], '\t', "%0.4f" % get_zscore(w, neg, H_neg)

"""
      T GT AG CT G C 		13 	8.1815 	5 	-0.7119
      T G AG CT GT C 		12 	7.4474 	6 	-0.3302
       T G A CT AG C 		8 	7.3425 	4 	0.3062
     GT CG AG CT G C 		17 	6.8985 	20 	1.6920
   GT CT AT A ACGT T 		26 	6.8194 	55 	5.2552
      T CG AG CT G C 		11 	6.7133 	3 	-1.4752
      T GT AG AC G C 		11 	6.7133 	6 	-0.3302
      T G AG C GT CT 		11 	6.7133 	4 	-1.0935
      T G AG CT CG C 		11 	6.7133 	2 	-1.8568
   T GT AT A ACGT CT 		25 	6.4523 	67 	7.5452
   A AT ACGT CT GT T 		25 	6.4523 	72 	8.4993
     GT T AG AC AT A 		16 	6.3794 	24 	2.7715
     GT CT AT A AG T 		16 	6.3794 	31 	4.6606
     GT CT AT A AC T 		16 	6.3794 	32 	4.9304
    T GT AG ACGT G C 		16 	6.3794 	9 	-1.2766
   T GT A AT ACGT CT 		24 	6.0853 	72 	8.4993
   AT ACGT A CT GT T 		24 	6.0853 	66 	7.3543
   T CT AT A ACGT GT 		24 	6.0853 	62 	6.5910
   CT AT ACGT T GT A 		24 	6.0853 	57 	5.6369
      T G AG AC GT C 		10 	5.9791 	4 	-1.0935
"""
