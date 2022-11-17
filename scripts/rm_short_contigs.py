import sys,os
dd = {}
len_filter = int(sys.argv[1]) #1000000
out_file = sys.argv[4] #"genome.masked.rm95.fa"
fout2 = open(out_file,'w')
with open(sys.argv[2],'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        if int(line_temp[1]) < len_filter:
            continue
        dd[line_temp[0]] = float(line_temp[1])
flag = 0
with open(sys.argv[3],'r') as fp:
    for line in fp:
        if line.startswith('>'):
            title = line.strip()[1:]
            if title not in dd:
                flag = 0
                continue
            fout2.write(line)
            flag = 1
        else:
            if flag == 1:
                fout2.write(line)
fout2.close()
       
os.system(f"samtools faidx {out_file}")













