import sys,os
dd = {}
ddd = {}
#with open("genome.masked.fa.fai",'r') as fp:
with open(sys.argv[1],'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        dd[line_temp[0]] = float(line_temp[1])
dict_fa = {}
#with open("genome.masked.fa",'r') as fp:
with open(sys.argv[2],'r') as fp:
    for line in fp:
        if line.startswith('>'):
            title = line.strip()[1:]
            ddd.setdefault(title,0)
            dict_fa.setdefault(title,[])
        else:
            chars = line.strip()
            masked = chars.count('a')+chars.count('c')+chars.count('g')+chars.count('t')+chars.count('N')+chars.count('n')
            ddd[title] += masked
            dict_fa[title].append(line)
out_file = sys.argv[3] #"genome.masked.rm95.fa"
#fout = open("masked_proportion.txt",'w')
fout = open(sys.argv[4],'w')
fout2 = open(out_file,'w')
for contig in dd:
    frac = ddd[contig]/dd[contig]
    fout.write(f'{contig}\t{int(dd[contig])}\t{round(frac,4)}\n')
    if dd[contig] * (1-frac) < 10000:
        continue
    fout2.write('>'+contig+'\n')
    fout2.write(''.join(dict_fa[contig]))
fout.close()
fout2.close()
       
os.system(f"samtools faidx {out_file}")













