library('SNFtool')

## First, set all the parameters:
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 10; 	# Number of Iterations, usually (10~20)

truelabel = c(matrix(1,53,1),matrix(2,54,1)); ##the ground truth of the simulated data;

W1<-as.matrix(read.table("../Scripts/3d_net.txt"))
#class(W1)
W2<-as.matrix(read.table("../Scripts/1d_net.txt"))


W = SNF(list(W1,W2), K, T)
write.table(W,file="../Results/1d_3d_fusionnet",sep='\t')

