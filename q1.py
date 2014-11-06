import re
from itertools import permutations
from sklearn.externals import joblib

promoter_pos = []
f = open('data/promoterPositive.fa.txt')
for line in f:
    if line[0] == '>':
        continue
    promoter_pos.append(line.strip())
f.close()

promoter_neg = []
f = open('data/promoterNegative.fa.txt')
for line in f:
    if line[0] == '>':
        continue
    promoter_neg.append(line.strip())
f.close()

L = ['A', 'C', 'G', 'T', '[AC]', '[AG]', '[AT]', '[CG]', '[CT]', '[GT]', '[ACGT]']
H_pos = {}
H_neg = {}
for i, p in enumerate(permutations(L, 6)):
    if i % 1000 == 0:
        print i
    patt = ''.join(p)
    H_pos[patt] = 0
    H_neg[patt] = 0
    for seq in promoter_pos:
        num_matches = len(re.findall(patt, seq))
        H_pos[patt] += num_matches
    for seq in promoter_neg:
        num_matches = len(re.findall(patt, seq))
        H_neg[patt] += num_matches

joblib.dump(H_pos, 'blobs/H_pos.pkl')
joblib.dump(H_neg, 'blobs/H_neg.pkl')
