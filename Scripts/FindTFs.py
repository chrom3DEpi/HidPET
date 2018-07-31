# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 09:56:12 2017

@author: WaRM
"""
import os
import commands
import argparse

##Requirements
##============
#Bedtools
#FIMO
#PASTAA
##============
                
#The Script runs a pipeline that uses the software FIMO(http://meme-suite.org/doc/fimo.html) and PASTAA(http://trap.molgen.mpg.de/PASTAA.htm) to scan for occurrences of a given motif in enhancer and promoter regions, and then find the confidence 3d TFs 

#Usage python FindTFs.py [options]
#      Options
#       -p   Path to the folder with HiPET data. Support both absolute
#       -i  [MANDATORY] Input file of the format "Enhancer.chr Enhancer.start Enhancer.end Promoter.chr Promoter.start Promoter.end¡° Interactions  
#       -r  [MANDATORY] Input file of your reference genome in fasta
#       -motifdbt  [MANDATORY] Input your motif db in transfac format
#       -motifdbm  [MANDATORY] Input your motif db in meme format
#       -pv  [OPTIONAL] Input your PASTAA P value set to choose your confidence TFs
#       -fv  [OPTIONAL] Input your FIMO P value set to choose your confidence TFs
#       -O   [OPTIONAL ]--output Your 3D TFs output file
##############

def getargs():
    ##Construct an ArgumentParser object for command-line arguments
    parser = argparse.ArgumentParser(description='HiPET FindTFs',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ##Path
    parser.add_argument('-p', '--path',default = '..',help = 'Path to the folder with HiPET data. Support both absolute')

    ##Input File
    parser.add_argument('-i', type=str, help='Input file of the format Enhancer.chr Enhancer.start Enhancer.end Promoter.chr Promoter.start Promoter.end Interactions',required=True)
    parser.add_argument('-r', type=str, help='Input file of your reference genome in fasta',required=True)
    parser.add_argument('-motifdbt', type=str, help='Input your motif db in transfac format ',required=True)
    parser.add_argument('-motifdbm', type=str, help='Input your motif db in meme format ',required=True)
    
    
    #P value
    parser.add_argument('-pv', type=float, help='Please input your PASTAA P value set to choose your confidence TFs',default=0.0001 )
    parser.add_argument('-fv', type=float, help='Please input your FIMO P value set to choose your confidence TFs',default=0.0001 )

#    Output File
    parser.add_argument('-O','--output', type=str, help='Your 3D TFs output file',default='tf_3d.list')

    return parser

def runfimo():
    parser = getargs()
    args = parser.parse_args()
    
    (status,fimo_outfil)=commands.getstatusoutput('fimo --oc fimo_output --parse-genomic-coord %s enhancer.fa' % args.motifdbm)
    os.system('sort fimo_output/fimo.txt |sort -n|uniq > fimo_unq.txt')
    
    
def fimo_tf():
    parser = getargs()
    args = parser.parse_args()
    
    with open('fimo_unq.txt','r') as infil:
        enhtf={}
        for line in infil:
            if line.startswith("#"):
                continue
            else:
                parse=line.rstrip().split()
                k=(parse[1],parse[2],parse[3])
                if float(parse[6])<float(args.fv):
                    if k in enhtf:
                        enhtf[k]=str(enhtf[k])+','+str('nrMotif'+parse[0])
                    else:
                        enhtf[k]='nrMotif'+parse[0]
    with open('fimo_unq_new.txt','w') as output:
        for i in enhtf:
            output.write(str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+str(enhtf[i])+'\n')
            
    
    os.system('intersectBed -a fimo_unq_new.txt -b enhancer.bed -wa -wb > fimo_enh_tf.txt')
    os.system('cut -f 4,5,6,7 fimo_enh_tf.txt > fimo_enh_tf_last.txt')
    with open('fimo_enh_tf_last.txt','r') as inf:
        global enh_tf
        enh_tf={}
        for line in inf:
            case=line.rstrip().split()
            enh_k=(case[1],case[2],case[3]) 
            if enh_k in enh_tf:
                enh_tf[enh_k]=enh_tf[enh_k]+','+str(case[0])
            else:
                enh_tf[enh_k]=str(case[0])


def pastaa_tf():
    parser = getargs()
    args = parser.parse_args()
    ##Compile
    os.system('chmod a+x %s/software/PASTAA/PSCM_to_PSEM' %args.path)
    os.system('chmod a+x %s/software/PASTAA/TRAP' %args.path)
    os.system('%s/software/PASTAA/PSCM_to_PSEM %s > %s/software/PASTAA/BRE_ENERGY_HOMO' %(args.path,args.motifdbt,args.path))
    os.system('%s/software/PASTAA/TRAP %s/software/PASTAA/BRE_ENERGY_HOMO promoter.fa  >  BRE_AFFY_HOMO' %(args.path,args.path))
    with open('BRE_AFFY_HOMO','r') as infil:
        global pasta_ans
        pasta_ans={}
        firstline=infil.readline()
        parse=firstline.rstrip().split()
        for line in infil:
            li=[]
            v=[]
            if line.startswith('ch'):
                case=line.rstrip().split()
                a=case[0].split(':')
                b=a[1].split('-')
                k=(str(a[0]),str(b[0]),str(b[1]))
                for i in case:
                    if i.startswith('ch'):
                        continue
                    else:
                        if float(i)<float(args.pv):
                            a=case.index(i)
                            li.append(a)
            if k in pasta_ans:
                continue
            else:                
                for j in li:
                    v.append(parse[j])
                pasta_ans[k]=v

def threeDtf():
    parser = getargs()
    args = parser.parse_args()
    os.chdir('../Results/')
    filename=args.output
    F=open(filename,'w')
    tf=[]
    
    with open(args.i,'r') as inf:
        loop=[]
        for line in inf:
            sourse=line.rstrip().split()
            e=(str(sourse[0]),str(sourse[1]),str(sourse[2]))
            p=(str(sourse[3]),str(sourse[4]),str(sourse[5]))
            loop.append([e,p])
    for i in loop:    
        if (i[0] in enh_tf.keys()) and (i[1] in pasta_ans.keys()):
            fimo_motif=enh_tf[i[0]].split(',')
            overlap=list(set(fimo_motif) & set(pasta_ans[i[1]]))
            for j in overlap:
                tf.append(j)
    tf={}.fromkeys(tf).keys()
    for t in tf:
        F.write(str(t)+'\n')     

    F.close()

            
def main():
    parser = getargs()
    args = parser.parse_args()
    
    
    if not os.path.isfile(args.i):
        parser.error("The input interaction file does not exist!")

    ##get the bed file from interactions file
    os.system('cut -f 1,2,3 %s > enhancer.bed' %args.i)
    os.system('cut -f 4,5,6 %s > promoter.bed' %args.i)
     
    ##get the fasta file from the bed file
    
    os.system('bedtools getfasta -fi %s -bed enhancer.bed -fo enhancer.fa' % args.r)
    os.system('bedtools getfasta -fi %s -bed promoter.bed -fo promoter.fa' % args.r)

    runfimo()
    fimo_tf()
    pastaa_tf()
    threeDtf()
    
    os.chdir('../Scripts/')
    os.remove('fimo_unq.txt')
    os.remove('fimo_unq_new.txt')
    os.remove('fimo_enh_tf.txt')
    os.remove('fimo_enh_tf_last.txt')
        

if __name__=='__main__':
    main()
    
#example
#python FindTFs.py -i ../FindTFs_example/interactions.txt -r ../Genome/hg19.fa -motifdbt ../FindTFs_example/motifdb_matrix.transfac -motifdbm ../FindTFs_example/motifdb_matrix.meme
