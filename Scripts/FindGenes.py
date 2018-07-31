# -*- coding: utf-8 -*-
"""
Created on Thu Sep 07 09:56:12 2017

@author: WaRM
"""
import os
import argparse

##Requirements
##============
#Bedtools
##============
                
#The Script can find the genes are regulated by maximal clique TFs

#Usage python FindGenes.py [options]
#      Options
#       -p   Path to the folder with HiPET data. Support both absolute
#       -i  [MANDATORY] input your species file of gene and promoter 
#       -c  [MANDATORY] Input Your cliques
#       -fv  [OPTIONAL] input your PASTAA P value set
#       -O   [OPTIONAL ]--output Your clique genes
##############



def getargs():
    parser = argparse.ArgumentParser(description='Parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i',type=str,help='Please input your species file of gene and promoter')
    parser.add_argument('-c',type=str,help='Your cliques',required=True)
    
    
    parser.add_argument('-fv', type=float, help='Please input your PASTAA P value set', required=True,default=0.0001 )

    parser.add_argument('-O','--output', type=str, help='Your clique',default='clique.gene')
    return parser




def Pro_genename():
    parser = getargs()
    args = parser.parse_args()
    
    global p_genename
    p_genename={}
    os.system('intersectBed -a promoter.bed -b %s -wa -wb > mypromoter.txt' % args.i)
    os.system('cut -f 1,2,3,7 mypromoter.txt > mypromoter_genename.txt')
    os.system('sort mypromoter_genename.txt |sort -n|uniq > mypromoter_genename_uniq.txt')
    os.system('rm mypromoter_genename.txt ')
    
                
    with open('mypromoter_genename_uniq_alone.txt','w') as output:
        with open('mypromoter_genename_uniq.txt','r') as infil:
            dic={}
            for line in infil:
                parse=line.rstrip().split()
                if (parse[0],parse[1],parse[2]) in dic:
                    v.append(parse[3])
                    dic[(parse[0],parse[1],parse[2])]=v
                else:
                    v=[]
                    v.append(parse[3])
                    dic[(parse[0],parse[1],parse[2])]=v
                    
            for i in dic:
                if len(dic[i])==1:
                    output.write(str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+str(dic[i][0])+'\n')

    with open('mypromoter_genename_uniq_alone.txt','r') as infil:
        for line in infil:
            parse=line.rstrip().split()
            if (parse[0],parse[1],parse[2]) in p_genename:
                p_genename[(parse[0],parse[1],parse[2])]=str(p_genename[(parse[0],parse[1],parse[2])])+'\t'+str(parse[3])
            else:
                p_genename[(parse[0],parse[1],parse[2])]=str(parse[3])

#
def P_ANS():
    parser = getargs()
    args = parser.parse_args()
    
    global pasta_ans
    pasta_ans={}
    with open('BRE_AFFY_HOMO','r') as infil:
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
                        if float(i)<float(args.fv):
                            a=case.index(i)
                            li.append(a)
            if k in pasta_ans:
                continue
            elif (k in p_genename):                
                for j in li:
                    v.append(parse[j])
                pasta_ans[p_genename[k]]=v

    #Read the 3d tf list from last step
    os.chdir('../Results/')
    with open('tf_3d.list','r') as inf:
        tf=[]
        for line in inf:
            parse=line.rstrip()
            tf.append(parse)
    os.chdir('../Scripts/')
         
    
    with open('genename_tf.txt','w') as output:
        for i in pasta_ans:
            output.write(str(i)+'\t')
            for j in pasta_ans[i]:
                if j in tf:
                    output.write(str(j))
            output.write('\n')
#            
def CliquetoGene(): 
    parser = getargs()
    args = parser.parse_args()
    
    with open('clique_gene.txt','w') as cliquegene:
        with open('genename_tf.txt','r') as infil:
            genametf={}
            for line in infil:
                parse=line.rstrip().split()
                new={}.fromkeys(parse).keys()
                genametf[parse[0]]=parse[1:len(new)]

        with open(args.c,'r') as inf:
            n=0
            for line in inf:
                n=n+1
                c='clique'+str(n)
                cliquegene.write(str(c)+'\n')
                c=line.rstrip().split()
                li=[]
                for i in genametf.keys():
                    if (len(list(set(c) & set(genametf[i]))) == len(c)) and (i not in li):
                        li.append(i)
                cliquegene.write('\t'.join(li)+'\n')
                
    os.system('rm genename_tf.txt')
    os.system('rm mypromoter.txt')
    os.system('rm mypromoter_genename_uniq.txt')
    os.system('rm mypromoter_genename_uniq_alone.txt')
#                    
def main():
    Pro_genename()
    P_ANS()
    CliquetoGene()
if __name__=='__main__':
    main()            
