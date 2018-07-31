#!/bin/bash
## HidPET
## HidPET (Hierarchical and dynamics analysis of TF cooperation with ChIA-PET and ChIP-Seq data) is a method to study the hierarchy and dynamics of ##TF cooperation  by integration of both ChIP-Seq and ChIA-PET datasets. 

## It contains 3 procedures:
## 1. find 3D TFs
## 2. Network fusion of 1D network and 3D TF network
## 3. clique detecting of TFs regulated genes.


echo ">>>>> HidPET "
echo ">>>>>  An Example for HidPET"

#basepath=$(cd `dirname $0`; pwd)
basepath="$(cd .. ; pwd)"
dir_FindTFs=$basepath/FindTFs_example
dir_Genome=$basepath/Genome
dir_NetworkFusion=$basepath/NetworkFusion_example
dir_FindGenes=$basepath/FindGenes_example

echo ">>>>> Step 1: TF_Finding "
TF_Finding="python FindTFs.py -i $dir_FindTFs/interactions.txt -r $dir_Genome/hg19.fa -motifdbt $dir_FindTFs/motifdb_matrix.transfac -motifdbm $dir_FindTFs/motifdb_matrix.meme"
$TF_Finding
echo $TF_Finding

echo ">>>>> Step 2: NetworkFusion "
NetworkFusion="python NetworkFusion.py -i $dir_NetworkFusion/tf_3d.list -I $dir_NetworkFusion/tf_1d.list -c $dir_NetworkFusion/human.ppi3d -C $dir_NetworkFusion/human.ppi1d"
$NetworkFusion
echo $NetworkFusion

echo ">>>>> Step 3: Clique_FindGenes "
Clique_FindGenes="python FindGenes.py -i $dir_FindGenes/promoter_genename -c $dir_FindGenes/maximal_clique -fv 0.0001 "
$Clique_FindGenes
echo $Clique_FindGenes




