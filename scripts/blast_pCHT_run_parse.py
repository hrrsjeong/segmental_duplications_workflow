import sys
fout = open(sys.argv[2],'w')
with open(sys.argv[1],'r') as fp:
    for line in fp:
        line_temp = line.strip().split('\t')
        if line.startswith("#"):continue
        if float(line_temp[-2]) > 1:continue
        fout.write(line_temp[1]+"\t"+str(min(int(line_temp[8]),int(line_temp[9])))+"\t"+str(max(int(line_temp[8]),int(line_temp[9])))+"\tpCHT\t"+line_temp[-2]+'\n')

