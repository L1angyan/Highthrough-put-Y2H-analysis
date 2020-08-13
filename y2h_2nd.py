#2020/07/29
#find PPI from 2nd NGS
#get no repeator and interactor like AB-BA

import sys
blastf1,blastf2,out = sys.argv[1],sys.argv[2],sys.argv[3]
readdict = {}
ppidict = {}
with open(blastf1) as f1:
    while 1:
        line1 = f1.readline().strip()
        if line1 == "":break
        if "#" in line1 or "_" not in line1:continue
        line1list = line1.split("\t")
        read,transcript1,align_length,mismatch = line1list[0],line1list[1],line1list[3],line1list[4]
        gene1 = transcript1[:transcript1.index("_")]
        if gene1 == 'T':
            readdict[read] = gene1
        elif int(align_length) >= 100 and int(mismatch) <= 5:
            readdict[read] = gene1
        while 1:
            line1 = f1.readline().strip()
            if "#" in line1:break

with open(blastf2) as f2:
    while 1:
        line2 = f2.readline().strip()
        if line2 == "":break
        if "#" in line2 or "_" not in line2:continue
        line2list = line2.split("\t")
        read,transcript2,align_length,mismatch = line2list[0],line2list[1],line2list[3],line2list[4]
        gene2 = transcript2[:transcript2.index("_")]
        indicator = 0
        if read in readdict and readdict[read] != gene2:
            if gene2 == "T":
                indicator = 1
            elif int(align_length) >= 100 and int(mismatch) <= 5:
                if readdict[read] == "T":
                    indicator = 1
                elif abs(int(readdict[read][-6:])-int(gene2[-6:])) > 5:
                    indicator = 1
        if indicator == 1:
            if readdict[read]+"\t"+gene2 in ppidict:
                ppidict[readdict[read]+"\t"+gene2]+=1
            elif gene2+"\t"+readdict[read] in ppidict:
                ppidict[gene2+"\t"+readdict[read]]+=1
            else:
                ppidict[readdict[read]+"\t"+gene2]=1
        while 1:
            line2  = f2.readline().strip()
            if "#" in line2:break

with open(out,'w') as op:
    for i,j in ppidict.items():
        op.write(i+"\t"+str(j)+"\n")
