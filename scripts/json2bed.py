import sys,os
import json

fout = open(sys.argv[2],'w')
f = open(sys.argv[1],'r')
data = json.load(f)
for repeat in data["Repeats"]:
    if "SequenceName" not in repeat:continue
    fout.write(f'{repeat["SequenceName"]}\t{repeat["Start"]}\t{int(repeat["Start"])+int(repeat["Length"])}\t{repeat["Consensus"]}\n')

f.close()
fout.close()
'''
"Repeats": [{"PassID": 0,
"SequenceName": "h1tg000008l",
"Start": 197081,
"Length": 1398,
"Period": 31,
"Score": 2063.160535,
"Log2 Pval": -377.192361,
"Substitutions": 46,
"Insertions": 0,
"Deletions": 2,
"Consensus": "CTGCCCCTTCCCCGACTGCCCCTTTCCATCT",
"Sequence": "CTGCCCCTTCCCCGACTGCCCCTTTCCATCTCCGCCCCTTCCCCGACTGCCCCTTTCCATCTCCGCCCCTTCCCCGACTGCCCCTTTCCATCTCCGCCCCTTCCCCGACTGCCCCTTTCCATCTCCGCCCCTTCCCCGACTGCCCCTTTCCATCTCCGCCCCTTCCCCGACTGCCCCTTTCCATCT
'''

'''
   470   14.6  2.3  2.4  chr10_RagTag                                2       650 (142677412) + (ACCCCC)n          Simple_repeat           1    663     (0)       1
  1704    1.4  2.1  2.3  chr10_RagTag                              651      2536 (142675526) + (CCCTAA)n          Simple_repeat           1   1882     (0)       2

chr18_RagTag    0       924     (GCAAG)n        924     +       Simple_repeat   unknown -1.0    1359391
chr3_RagTag     0       97      MER21A  97      +       LTR     ERVL    -1.0    2377419
'''
