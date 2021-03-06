HidPET
Introduction
HidPET (Hierarchical and dynamics analysis of TF cooperation with ChIA-PET and ChIP-Seq data) is a method to study the hierarchy and dynamics of TF cooperation by integration of both ChIP-Seq and ChIA-PET datasets. This package contains the  scripts and examples to find 3D TFs, Network fusion of 1D network and 3D TF network, clique detecting of TFs regulated genes.
Three procedures are given below which are TFs Finding, Network Fusion and Clique detecting of TFs regulated genes. 

How to run it ?
HidPET is developed in python and can be downloaded from https://github.com/chrom3DEpi/HidPET. This repository contains scripts,examples and required software for HidPET.

  Scripts
		FindTFs.py
		NetworkFusion.py
		FindGenes.py
  Software
		PASTAA // Find the TFs in promoter region
		net_fusion.R // network fusion pipeline in R
		SNFtool-master.zip // Package of network fusion

When you run the program, please change your working directory to HidPET/Scripts/. 

In addition, HidPET requires the following dependencies: 
##Bedtools (https://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz)
#FIMO (http://meme-suite.org/doc/fimo.html)
#PASTAA (http://trap.molgen.mpg.de/PASTAA.htm)
#SNFtool (https://cran.r-project.org/src/contrib/SNFtool_2.2.1.tar.gz)
#R (http://www.r-project.org/)
Please install the above dependencies and then read the README file for more information!

Required data
To run HidPET, the following data should be prepared:
(i)the reference genome sequence data, as an example hg19 (ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz)
(ii)chromatin interaction location which is BEDPE format (http://bedtools.readthedocs.io/en/latest/content/general-usage.html), 
(iii)PWM for TF motif in meme and transfac format.

Execution
HidPET is an easy-to-use pipeline and it is run in the following command. 
sh run_clique.sh


Examples:
An example script for each procedure is given as the following:
1. TFs Finding
python FindTFs.py -i ../FindTFs_example/interactions.txt -r ../Genome/hg19.fa -motifdbt ../FindTFs_example/motifdb_matrix.transfac -motifdbm ../FindTFs_example/motifdb_matrix.meme

2. Network Fusion
python NetworkFusion.py -i ../NetworkFusion_example/tf_3d.list -I ../NetworkFusion_example/tf_1d.list -c ../NetworkFusion_example/human.ppi3d -C ../NetworkFusion_example/human.ppi1d

3. Clique detecting of TFs regulated genes
python FindGenes.py -i ../FindGenes_example/promoter_genename -c ../FindGenes_example/maximal_clique -fv 0.0001 


Contacct us
If you have any questions or suggestions, please send an email to Ruimin Wang(rmwang@webmail.hzau.edu.cn).

---END---





