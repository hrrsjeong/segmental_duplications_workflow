import sys,os
fout = open(sys.argv[2],'w')
with open(sys.argv[1],'r') as fp:
    cnt = 0
    fp.readline()
    fp.readline()
    fp.readline()
    for line in fp:
        line_temp = [x.strip() for x in line.strip().split() if x.strip() != '']
        #print (line_temp)
        #cnt +=1
        #if cnt == 10:break
        chrom,start,end = line_temp[4],str(int(line_temp[5])-1),line_temp[6]
        rtype,rlen,strand = line_temp[9],int(end)-int(start),"+"
        rr = line_temp[10].split('/')
        if len(rr) == 2:
            rfam,rclass = rr
        elif len(rr) == 1:
            rfam,rclass = rr[0],"Unknown"
        else:
            print ("repeat class error")
            break
        
        fake_int1,fake_int2 = "-1","1000"
        fout.write(f"{chrom}\t{start}\t{end}\t{rtype}\t{rlen}\t{strand}\t{rfam}\t{rclass}\t{fake_int1}\t{fake_int2}\n")
fout.close()

'''
   470   14.6  2.3  2.4  chr10_RagTag                                2       650 (142677412) + (ACCCCC)n          Simple_repeat           1    663     (0)       1
  1704    1.4  2.1  2.3  chr10_RagTag                              651      2536 (142675526) + (CCCTAA)n          Simple_repeat           1   1882     (0)       2

chr18_RagTag    0       924     (GCAAG)n        924     +       Simple_repeat   unknown -1.0    1359391
chr3_RagTag     0       97      MER21A  97      +       LTR     ERVL    -1.0    2377419
'''
