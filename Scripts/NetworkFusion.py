# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:12:48 2017

@author: WaRM
"""
import os
import argparse

##Requirements
##============
#R
#SNFtool(https://cran.r-project.org/src/contrib/SNFtool_2.2.1.tar.gz)
##============
                
#The Script runs a pipeline that uses the R packege "SNFtool" to fusion the 1d network and 3d network

#Usage python NetworkFusion.py [options]
#      Options
#       -p   Path to the folder with HiPET data. Support both absolute
#       -i  [MANDATORY] Input file of your 3d TFs list  
#       -I  [MANDATORY] Input file of your 1d TFs list
#       -c [MANDATORY] Input file of your TFs interactions file in 3d
#       -C  [MANDATORY] Input file of your TFs interactions file in 1d
#       -o   [OPTIONAL ]--output Your fusion network file
##############

def getargs():
    ##Construct an ArgumentParser object for command-line arguments
    parser = argparse.ArgumentParser(description='HiPET NetworkFusion',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    ##Path
    parser.add_argument('-p', '--path',default = '..',help = 'Path to the folder with HiPET data. Support both absolute')
    
    ##Input File
    parser.add_argument('-i', type=str, help='Input file of your 3d TFs list',required=True)
    parser.add_argument('-I', type=str, help='Input file of your 1d TFs list',required=True)
    parser.add_argument('-c', type=str, help='Input file of your TFs interactions file in 3d',required=True)
    parser.add_argument('-C', type=str, help='Input file of your TFs interactions file in d',required=True)
    parser.add_argument('-o', type=str, help='Your fusion network file',default='1d_3d_fusionnet')
    return parser


def TF_NET():
    parser = getargs()
    args = parser.parse_args()
    with open('3d_tf_tf_cor.txt','w') as output:
        with open(args.i,'r') as inf:
            li=[]
            for line in inf:
                 parse=line.rstrip().split()
                 li.append(parse[0])
        with open(args.c,'r') as infil:
            for line in infil:
                case=line.rstrip().split()
                if (case[0] in li) and (case[1] in li) and int(case[2])>400:
                    output.write(line)

def TFoverlap():
    parser = getargs()
    args = parser.parse_args()
    with open('tf_overlap.txt','w') as output:
        with open(args.i,'r') as inf:
            li1=[]
            for line in inf:
                parse=line.rstrip().split()
                li1.append(parse[0])
    
        with open(args.I,'r') as infil:
            li2=[]
            for line in infil:
                case=line.rstrip().split()
                li2.append(case[0])
        overlap=list(set(li1) & set(li2))
        for i in overlap:
            output.write(str(i)+'\n')
            
def TF_t_net():
    parser = getargs()
    args = parser.parse_args()
    d3=open('3d_net.txt','w')
         
    with open('tf_overlap.txt','r') as inf:
        li=[]
        for line in inf:
            case=line.rstrip().split()
            li.append(case[0])
            
        d3.write('\t'+'\t'.join(li)+'\n')

        table=[[0 for i in range(len(li))] for j in range(len(li))]
    
    
    with open(args.c,'r') as infil:
        for line in infil:
            parse=line.rstrip().split()
            if (parse[0] in li) and (parse[1] in li):
                table[li.index(parse[0])][li.index(parse[1])]=parse[2]
        
    for i in range(len(li)):
        d3.write(str(li[i]))
        for j in range(len(li)):
            if i==j:
                d3.write('\t'+'1000')
            else:
                d3.write('\t'+str(table[i][j]))
        d3.write('\n')
    
def TF_o_net():
    parser = getargs()
    args = parser.parse_args()
    
    d1=open('1d_net.txt','w')      
    with open('tf_overlap.txt','r') as inf:
        li=[]
        for line in inf:
            case=line.rstrip().split()
            li.append(case[0])

        d1.write('\t'+'\t'.join(li)+'\n')
        table=[[0 for i in range(len(li))] for j in range(len(li))]
           
    with open(args.C,'r') as infil:
        for line in infil:
            parse=line.rstrip().split()
            if (parse[0] in li) and (parse[1] in li):
                table[li.index(parse[0])][li.index(parse[1])]=parse[2]
           
    for i in range(len(li)):
        d1.write(str(li[i]))
        for j in range(len(li)):
            d1.write('\t'+str(table[i][j]))
        d1.write('\n')

def NET_fusion():
    parser = getargs()
    args = parser.parse_args()
    
    os.system('Rscript ../software/net_fusion.R')
    os.system('sed -i "1d" ../Results/%s' %args.o)
    os.system('perl -pi -e "s/.+?\t//" ../Results/%s' %args.o)
    
    
def main():
    getargs()
    TF_NET()
    TFoverlap()
    TF_t_net()
    TF_o_net()
    NET_fusion()
    
    os.remove('3d_tf_tf_cor.txt')
    os.remove('1d_net.txt')
    os.remove('3d_net.txt')
    os.remove('tf_overlap.txt')
    
if __name__=='__main__':
    main()
    
#example
#python NetworkFusion.py -i ../NetworkFusion_example/tf_3d.list -I ../NetworkFusion_example/tf_1d.list -c ../NetworkFusion_example/human.ppi3d -C ../NetworkFusion_example/human.ppi1d
