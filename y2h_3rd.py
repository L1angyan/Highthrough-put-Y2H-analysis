# -*- coding: utf-8 -*-
#python
#find pais of ad_gene and bd_gene by squence
#Authors:Liangyan in HuaZhong Agriculture University/WuHan/Hubei/China
#2019/08/08

import os
#os.system("linux command") execute command
#a = os.popen("linux command") assignment
import sys #recive outside parameter
import numpy as np

fafile = sys.argv[2]
outputfile = sys.argv[3]
DETAIL = outputfile + '_detail'
os.system(" ".join(("formatdb -i" , fafile , "-n" , fafile , "-o T"))) #将fafile格式化，作为比对序列数据库

def judgeOverlap(startlist,endlist,blastlineline):       
    blastlinelinecontentlist = blastlineline.split("\t")
    i=0
    res=0
    while(i<len(startlist)):#$i<10 而$res可能等于0，或者是$i+1
        if(endlist[i]!=0):
            if((startlist[i]>=int(blastlinelinecontentlist[6]) and endlist[i]<=int(blastlinelinecontentlist[7])) or (startlist[i]<=int(blastlinelinecontentlist[6]) and endlist[i]>=int(blastlinelinecontentlist[7])) or (startlist[i]>=int(blastlinelinecontentlist[6]) and startlist[i]<=int(blastlinelinecontentlist[7]) and (int(blastlinelinecontentlist[7])-startlist[i])/(endlist[i]-startlist[i])>0.3) or (startlist[i]>=int(blastlinelinecontentlist[6]) and startlist[i]<=int(blastlinelinecontentlist[7]) and (int(blastlinelinecontentlist[7])-startlist[i])/(int(blastlinelinecontentlist[7])-int(blastlinelinecontentlist[6]))>0.3) or (endlist[i]>=int(blastlinelinecontentlist[6]) and endlist[i]<=int(blastlinelinecontentlist[7]) and (endlist[i]-int(blastlinelinecontentlist[6]))/(endlist[i]-startlist[i])>0.3) or (endlist[i]>=int(blastlinelinecontentlist[6]) and endlist[i]<=int(blastlinelinecontentlist[7]) and (endlist[i]-int(blastlinelinecontentlist[6]))/(int(blastlinelinecontentlist[7])-int(blastlinelinecontentlist[6]))>0.3)):
            #检查是否有overlap的情况，有overlap再看其是否比对得比已知的更好。
                res=i+1
                break #一旦满足这个条件那么就找到了这个索引值就直接跳出这个while循环。
            
        
        i+=1
    
    #print "$res\n";
    return res

def parseNeedle(needlefile):
    firststr,secondstr = "",""
    score,start,end = 0,0,0
    IN = open(needlefile)
    tstr = 1
    while(tstr):
        tstr = IN.readline().strip()
        if ("# Score:" in tstr):
            score = float(tstr[tstr.index(":")+2:]) #读取score的值
        elif(">" in tstr) and (firststr == ""):
            while(tstr):
                tstr = IN.readline().strip()
                if(">" in tstr): break
                firststr += tstr.strip()
            
            while(tstr):
                tstr = IN.readline().strip()
                if("#" in tstr): break
                secondstr += tstr.strip()

#firststr是read序列
#secondstr序列是比对结果，未配对的序列用字符串“-”代替。

    i=0
    while(i<len(secondstr)):
        if(secondstr[i] != "-"):
            start += 1
            break
        else:
            if(firststr[i] != "-"):
                start += 1

        i+=1
    
    i = len(secondstr)-1
    end = len(secondstr)
    
    while(i >= 0):
        if(secondstr[i] != "-"):    
            break
        else:
            if(firststr[i] != "-"):
                end-=1
        i-=1
    
    IN.close()
    return [score,start,end]

def getLibInfo(targetfile,ADfile,ADRCfile,BDfile,BDRCfile,attLfile,attLRCfile):
    ADCMDstr = " ".join(("needle -gapopen 10 -gapextend 0.5",targetfile,ADfile,"-outfile ad -auto -adesshow3 Y -aformat3 markx3"))
    ADRCCMDstr = " ".join(("needle -gapopen 10 -gapextend 0.5",targetfile,ADRCfile,"-outfile adr -auto -adesshow3 Y -aformat3 markx3"))
    BDCMDstr = " ".join(("needle -gapopen 10 -gapextend 0.5",targetfile,BDfile,"-outfile bd -auto -adesshow3 Y -aformat3 markx3"))
    BDRCCMDstr = " ".join(("needle -gapopen 10 -gapextend 0.5",targetfile,BDRCfile,"-outfile bdr -auto -adesshow3 Y -aformat3 markx3"))
    attLCMDstr = " ".join(("needle -gapopen 10 -gapextend 0.5",targetfile,attLfile,"-outfile attL -auto -adesshow3 Y -aformat3 markx3"))
    attLRCCMDstr = " ".join(("needle -gapopen 10 -gapextend 0.5",targetfile,attLRCfile,"-outfile attLr -auto -adesshow3 Y -aformat3 markx3"))
    os.system(ADCMDstr)
    os.system(ADRCCMDstr)
    os.system(BDCMDstr)
    os.system(BDRCCMDstr)
    os.system(attLCMDstr)
    os.system(attLRCCMDstr)
    os.system("grep -v '^$' ad > ad1")
    os.system("grep -v '^$' adr > adr1")
    os.system("grep -v '^$' bd > bd1")
    os.system("grep -v '^$' bdr > bdr1")
    os.system("grep -v '^$' attL > attL1")
    os.system("grep -v '^$' attLr > attLr1")
    ADarr=parseNeedle("ad1")
    ADRCarr=parseNeedle("adr1")
    BDarr=parseNeedle("bd1")
    BDRCarr=parseNeedle("bdr1")
    attLarr=parseNeedle("attL1")
    attLRCarr=parseNeedle("attLr1")
    os.system("rm ad ad1")
    os.system("rm adr adr1")
    os.system("rm bd bd1")
    os.system("rm bdr bdr1")
    os.system("rm attL attL1")
    os.system("rm attLr attLr1")
    
    ADscore,ADstart,ADend = 0,0,0
    if(ADarr[0]<30 and ADRCarr[0]<30):
        pass
    elif(ADarr[0]<ADRCarr[0]):
        ADscore,ADstart,ADend = ADRCarr[0],ADRCarr[1],ADRCarr[2]
    else:
        ADscore,ADstart,ADend = ADarr[0],ADarr[1],ADarr[2]
    
    BDscore,BDstart,BDend = 0,0,0
    if(BDarr[0]<30 and BDRCarr[0]<30):
        pass
    elif(BDarr[0]<BDRCarr[0]):
        BDscore,BDstart,BDend = BDRCarr[0],BDRCarr[1],BDRCarr[2]
    else:
        BDscore,BDstart,BDend = BDarr[0],BDarr[1],BDarr[2]
    
    attLscore,attLstart,attLend = 0,0,0
    if(attLarr[0]<30 and attLRCarr[0]<30):
        pass
    elif(attLarr[0]<attLRCarr[0]):
        attLscore,attLstart,attLend = attLRCarr[0],attLRCarr[1],attLRCarr[2]
    else:
        attLscore,attLstart,attLend = attLarr[0],attLarr[1],attLarr[2]
    return [ADscore,ADstart,ADend,BDscore,BDstart,BDend,attLscore,attLstart,attLend]

def checkVectorStruct(positionlist,start1,end1,start2,end2):
    resflag=0
    
    if(positionlist[0]==0 or positionlist[3]==0 or positionlist[6]==0):
        pass
    else:
        Amid=(positionlist[1]+positionlist[2])/2
        Bmid=(positionlist[4]+positionlist[5])/2
        Mmid=(positionlist[7]+positionlist[8])/2
        Cmid=(end1+start1)/2
        Dmid=(end2+start2)/2
        if(Amid<=Cmid and Cmid<=Mmid and Mmid<=Dmid and Dmid<=Bmid):
            resflag=1 #ACMDB
        elif(Amid<=Dmid and Dmid<=Mmid and Mmid<=Cmid and Cmid<=Bmid):
            resflag=2 #ADMCB
        elif(Bmid<=Cmid and Cmid<=Mmid and Mmid<=Dmid and Dmid<=Amid):
            resflag=3 #BCMDA
        elif(Bmid<=Dmid and Dmid<=Mmid and Mmid<=Cmid and Cmid<=Amid):
            resflag=4 #BDMCA

    return resflag

blast_txt = open(sys.argv[1])
detail = open(DETAIL,'a')
detail.write('T3rdNGSRead\tAD\tBD\tOrderType\tADPrimerscore\tADPrimerstart\tADPrimerend\tBDPrimerscore\tBDPrimerstart\tBDPrimerend\tattLRecscore\tattLRecstart\tattLRecend\tR1start\tR1end\tG1start\tG1end\tR2start\tR2end\tG2start\tG2end\n')
blastline = 1 
ppidictionary = {}
ADfile,ADRCfile,BDfile,BDRCfile,attLfile,attLRCfile = "AD.fa","AD_revcom.fa","BD.fa","BD_revcom.fa","attL.fa","attL_revcom.fa"

while blastline:  
    blastline = blast_txt.readline().strip("\n")
    if "#" in blastline:continue
    if blastline == "":break
    blastlinecontentlist = blastline.split("\t") #blast结果每一行按列分开，写入list，方便取出元素
    readname = blastlinecontentlist[0]
    
    #ok 选reads比对上的可能是gene1 ，gene2的候选基因，通过judgeoverlap
    genelist = ["","","","","","","","","",""]
    startlist = [0,0,0,0,0,0,0,0,0,0]
    endlist = [0,0,0,0,0,0,0,0,0,0]
    Genestartlist = [0,0,0,0,0,0,0,0,0,0]
    Geneendlist = [0,0,0,0,0,0,0,0,0,0]
    scorelist = [0,0,0,0,0,0,0,0,0,0] #分别是基因ID，配对序列位于reads上的起始、终止，配对序列位于gene上的起始、终止，比对得分的列表
    genelist[0]=blastlinecontentlist[1]
    startlist[0]=int(blastlinecontentlist[6])
    endlist[0]=int(blastlinecontentlist[7])
    Genestartlist[0]=int(blastlinecontentlist[8])
    Geneendlist[0]=int(blastlinecontentlist[9])
    scorelist[0]=float(blastlinecontentlist[11])
    listi = 0 #等下judgeoverlap之后会用到
    blastlineline = 1
    while blastlineline: #承接上面继续逐行遍历没有#的小单元
        blastlineline = blast_txt.readline().strip("\n")
        if "#" in blastlineline: break #遇到#就跳出循环，保证处理的所有reads都属于一个单元内
        blastlinelinecontentlist = blastlineline.split("\t") #按列分隔每一行，方便取出元素 #ok数据读入完毕，最多选择10个每个reads匹配上的基因，原则:有overlap的话，留下overlap最好的那个配对，并存入列表。没有overlap的话，直接存入
        index = judgeOverlap(startlist,endlist,blastlineline)
        if index != 0:
            if scorelist[index-1] < float(blastlinelinecontentlist[11]): #不仅仅overlap要匹配上而且要score要比原先要高
                scorelist[index-1] = float(blastlinelinecontentlist[11])
                startlist[index-1] = int(blastlinelinecontentlist[6])
                endlist[index-1] = int(blastlinelinecontentlist[7])
                genelist[index-1] = blastlinelinecontentlist[1]
                Genestartlist[index-1] = int(blastlinelinecontentlist[8])
                Geneendlist[index-1] = int(blastlinelinecontentlist[9])
        else:#read比对到的基因与其他基因没有overlap，将其直接放入候选的基因列表中,并将候选的基因的关键信息存入列表
            listi += 1
            if listi == 10 :break
            genelist[listi] = blastlinelinecontentlist[1]
            startlist[listi] = int(blastlinelinecontentlist[6])
            endlist[listi] = int(blastlinelinecontentlist[7])
            scorelist[listi] = float(blastlinelinecontentlist[11])
            Genestartlist[listi] = int(blastlinelinecontentlist[8])
            Geneendlist[listi] = int(blastlinelinecontentlist[9])
    #得到长度为10的几个列表，其中包含基因ID，匹配序列在read上的起点终点，在基因上的起点终点，匹配序列得分（bit score）
    #现在开始对候选基因列表进行排序，选出候选的与ad连接的gene1，和与bd连接的gene2
    sortedindex = np.argsort(np.array(scorelist)).tolist()[::-1] #将列表按照score的值从大到小排序并且返回其索引值
    if (genelist[sortedindex[0]] != "") and (genelist[sortedindex[1]] != ""):
        readfile = "readfile.fa"
        os.system(" ".join(("fastacmd -d", fafile ,"-s",readname, ">" , readfile))) #fafile比对建库，比对read名以fa格式输出read序列
        positionlist = getLibInfo(readfile,ADfile,ADRCfile,BDfile,BDRCfile,attLfile,attLRCfile)
        os.system(" ".join(("rm",readfile)))
        vectType = checkVectorStruct(positionlist,startlist[sortedindex[0]],endlist[sortedindex[0]],startlist[sortedindex[1]],endlist[sortedindex[1]])
        genea,geneb = genelist[sortedindex[0]],genelist[sortedindex[1]]
        detail.write(readname+'\t'+genea+'\t'+geneb+'\t'+str(vectType)+'\t'+'\t'.join([str(i) for i in positionlist])+'\t'+str(startlist[sortedindex[0]])+'\t'+str(endlist[sortedindex[0]])+'\t'+str(Genestartlist[sortedindex[0]])+'\t'+str(Geneendlist[sortedindex[0]])+'\t'+str(startlist[sortedindex[1]])+'\t'+str(endlist[sortedindex[1]])+'\t'+str(Genestartlist[sortedindex[1]])+'\t'+str(Geneendlist[sortedindex[1]])+'\n')
        if(vectType):
            tkey = ""
            if(vectType==2 or vectType==3):
                tkey = genea +"\t"+ geneb
                
            else:
                tkey=geneb +"\t"+ genea
            
            if tkey in ppidictionary.keys():
                ppidictionary[tkey] += 1
                #如果该配对存在，hits数加一
            else:
                ppidictionary[tkey] = 1
#最后得到了一个字典，其键是两个配对的基因，其值为hits数
blast_txt.close()
detail.close()
f = open(outputfile,"a")
for i,j in ppidictionary.items():
    f.write(i+"\t"+str(j)+"\n")
f.close()
