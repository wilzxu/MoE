rm(list = ls())
library(ggplot2)
library(mixtools)
path = getwd()
setwd(path)
source("./multiplot.R")
### plot NOdata in mixtools pkg

###############NOdata#####################################
data("NOdata")
H_NO = read.table('output_H_NO.txt')
dataNO = data.frame("x" = NOdata[,1], "y" = NOdata[,2], "h1" = H_NO[,1], "h2" = H_NO[,2])
ggplot(NOdata, aes(x = NO, y = Equivalence )) + geom_point() +theme_minimal()
true_pars_NO = read.table('output_parameters_NO.txt')
true_pars_NO = c(t(as.matrix(true_pars_NO)))
est_pars_NO = read.table('output_log_NO.txt') 


p1 = qplot(1:nrow(est_pars_NO), est_pars_NO[,1], data = est_pars_NO, xlab = "", ylab = "expert1", main = "beta0") 
 # geom_hline(yintercept = true_pars_NO[1], color="red")
p2 = qplot(1:nrow(est_pars_NO), est_pars_NO[,2], data = est_pars_NO, xlab = "", ylab = "expert2") 
  #geom_hline(yintercept = true_pars_NO[2], color="red")
p3 = qplot(1:nrow(est_pars_NO), est_pars_NO[,3], data = est_pars_NO, xlab = "", ylab = "", main = "beta1") 
  #geom_hline(yintercept = true_pars_NO[3], color="red")
p4 = qplot(1:nrow(est_pars_NO), est_pars_NO[,4], data = est_pars_NO, xlab = "", ylab = "") 
  #geom_hline(yintercept = true_pars_NO[4], color="red")
p5 = qplot(1:nrow(est_pars_NO), est_pars_NO[,5], data = est_pars_NO, xlab = "", ylab = "", main = "theta0") 
  #geom_hline(yintercept = true_pars_NO[5], color="red")
p6 = qplot(1:nrow(est_pars_NO), est_pars_NO[,6], data = est_pars_NO, xlab = "", ylab = "") 
  #geom_hline(yintercept = true_pars_NO[6], color="red")
p7 = qplot(1:nrow(est_pars_NO), est_pars_NO[,7], data = est_pars_NO, xlab = "", ylab = "", main = "theta1") 
  #geom_hline(yintercept = true_pars_NO[7], color="red")
p8 = qplot(1:nrow(est_pars_NO), est_pars_NO[,8], data = est_pars_NO, xlab = "", ylab = "") 
  #geom_hline(yintercept = true_pars_NO[8], color="red")
p9 = qplot(1:nrow(est_pars_NO), est_pars_NO[,9], data = est_pars_NO, xlab = "", ylab = "", main = "sigma") 
  #geom_hline(yintercept = true_pars_NO[9], color="red")
p10 = qplot(1:nrow(est_pars_NO), est_pars_NO[,10], data = est_pars_NO, xlab = "", ylab = "") 
  #geom_hline(yintercept = true_pars_NO[10], color="red")

multiplot(p1, p2, p3, p4, p5,p6,p7,p8,p9,p10,cols=5)

h_plot = ggplot(dataNO, aes(x, y = value, color = variable)) + 
  geom_point(aes(y = h1, col = "h1")) + 
  geom_point(aes(y = h2, col = "h2")) +
  ggtitle("H")


likelihood_plot = qplot(1:nrow(est_pars_NO), est_pars_NO[,11], data = est_pars_NO, xlab = "", ylab = "", main = "loglikelihood")

multiplot(h_plot, likelihood_plot,cols=2)


#########j2################################
x2 = read.table('X_j2.txt')
y2 = read.table('Y_j2.txt')
H_2 = read.table('output_H_vN_j2.txt')
H_2_qN = read.table("output_H_qN_j2.txt")
H_2_gd = read.table('output_H_gd_j2.txt')
data2 = data.frame("int" = x2[,1], "x" = x2[,2], "y" = y2[,1], "id" = y2[,2], "h1" = H_2[,1], "h2" = H_2[,2])
ggplot(data2, aes(x = x, y = y )) + geom_point()



# 
# true_pars_2 = read.table('output_parameters_qNR_j2.txt')
# true_pars_2 = c(t(as.matrix(true_pars_2)))
# true_pars_2_gd = read.table('output_parameters_gd_j2.txt')
# true_pars_2_gd = c(t(as.matrix(true_pars_2_gd)))

maxiter = 1000
est_pars_2 = read.table('output_log_vN_j2.txt')
est_pars_2 = est_pars_2[1:maxiter,][-(1:10), ]

est_pars_2_qN = read.table('output_log_qN_j2.txt')
est_pars_2_qN = est_pars_2_qN[1:maxiter,][-(1:10), ]

est_pars_2_gd = read.table('output_log_gD_j2.txt')
est_pars_2_gd = est_pars_2_gd[1:maxiter,][-(1:10), ]




p1 = ggplot(data = est_pars_2) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,1], alpha = 0.1), data = est_pars_2_gd,  color = "grey") +
    geom_point(aes(1:nrow(est_pars_2), est_pars_2[,1]), alpha = 0.1, color = "darkred") +    
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,1], alpha = 0.1), data = est_pars_2_qN,  color = "darkblue") +
  geom_hline(yintercept = 2.5, color="red") +
   ggtitle("beta0") + xlab("iter") + ylab("expert1") + theme(legend.position="none") 

p2 = ggplot(data = est_pars_2) + theme_minimal() +   
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,2], alpha = 0.1), data = est_pars_2_gd, color = "grey") + 
  geom_point(aes(1:nrow(est_pars_2), est_pars_2[,2], alpha = 0.1), data = est_pars_2,  color = "darkred") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,2]), alpha = 0.1, data = est_pars_2_qN, color = "darkblue") +
  geom_hline(yintercept = 0.4, color="red") + 
  ylab("expert2") +
  ggtitle("beta0") + xlab("iter") + theme(legend.position="none")


p3 = ggplot(data = est_pars_2) + theme_minimal() +   
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,3], alpha = 0.1), data = est_pars_2_gd, color = "grey") + 
  geom_point(aes(1:nrow(est_pars_2), est_pars_2[,3], alpha = 0.1), data = est_pars_2,  color = "darkred") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,3]), alpha = 0.1, data = est_pars_2_qN, color = "darkblue") +
  geom_hline(yintercept = -0.5, color="red") + 
  ggtitle("beta1") + xlab("iter") + ylab("") + theme(legend.position="none")


p4 = ggplot(data = est_pars_2) + theme_minimal() +   
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,4], alpha = 0.1), data = est_pars_2_gd, color = "grey") + 
  geom_point(aes(1:nrow(est_pars_2), est_pars_2[,4], alpha = 0.1), data = est_pars_2,  color = "darkred") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,4]), alpha = 0.1,data = est_pars_2_qN, color = "darkblue") +
  geom_hline(yintercept = 0.8, color="red") +
  xlab("iter") + ylab("") + theme(legend.position="none")


p5 = ggplot(data = est_pars_2) + theme_minimal() +   
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,5], alpha = 0.1), data = est_pars_2_gd, color = "grey") + 
  geom_point(aes(1:nrow(est_pars_2), est_pars_2[,5], alpha = 0.1), data = est_pars_2,  color = "darkred") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,5]), alpha = 0.1, data = est_pars_2_qN, color = "darkblue") +
  ggtitle("theta0") +
 xlab("iter") + ylab("") + theme(legend.position="none")


p6 = ggplot(data = est_pars_2) + theme_minimal() +   
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,6], alpha = 0.1), data = est_pars_2_gd, color = "grey") + 
  geom_point(aes(1:nrow(est_pars_2), est_pars_2[,6], alpha = 0.1), data = est_pars_2,  color = "darkred") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,6]), alpha = 0.5, data = est_pars_2_qN, color = "darkblue") +
  xlab("iter") + ylab("") + theme(legend.position="none")


p7 = ggplot(data = est_pars_2) + theme_minimal() +   
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,7], alpha = 0.1), data = est_pars_2_gd, color = "grey") + 
  geom_point(aes(1:nrow(est_pars_2), est_pars_2[,7], alpha = 0.1), data = est_pars_2,  color = "darkred") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,7]), alpha = 0.1, data = est_pars_2_qN,  color = "darkblue") +
  ggtitle("theta1") +
 xlab("iter") + ylab("") + theme(legend.position="none")



p8 = ggplot(data = est_pars_2) + theme_minimal() +   
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,8], alpha = 0.1), data = est_pars_2_gd, color = "grey") + 
  geom_point(aes(1:nrow(est_pars_2), est_pars_2[,8], alpha = 0.1), data = est_pars_2,  color = "darkred") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,8]), alpha = 0.1,  data = est_pars_2_qN,  color = "darkblue") +
  xlab("iter") + ylab("") + theme(legend.position="none")


p9 = ggplot(data = est_pars_2) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,9], alpha = 0.1), data = est_pars_2_gd, color = "grey") +  
  geom_point(aes(1:nrow(est_pars_2), est_pars_2[,9], alpha = 0.1), data = est_pars_2,  color = "darkred") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,9]), alpha = 0.1,  data = est_pars_2_qN,  color = "darkblue") +
   geom_hline(yintercept = 0.09, color="red") + 
   xlab("iter") + ylab("") + ggtitle("sigma2") +theme(legend.position="none")

p10 = ggplot(data = est_pars_2) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,10], alpha = 0.1), data = est_pars_2_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2[,10], alpha = 0.1), data = est_pars_2,  color = "darkred") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,10]), alpha = 0.1, data = est_pars_2_qN,  color = "darkblue") +
  geom_hline(yintercept = 0.09, color="red") + 
  xlab("iter") + ylab("") + theme(legend.position="none")

multiplot(p1, p2, p3 , p4, p5, p6, p7, p8, p9 , p10,cols=5)


################
h_plot = ggplot(data2, aes(x, y = value, color = variable)) + 
  geom_point(aes(y = h1, col = "h1")) + 
  geom_point(aes(y = h2, col = "h2")) +
  ggtitle("H") +theme_minimal()


likelihood_plot = ggplot(data = est_pars_2) + 
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_gd[,11]),  alpha = 0.2, data = est_pars_2_gd, color ="grey") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2[,11]), data = est_pars_2, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_2), est_pars_2_qN[,11]), data = est_pars_2_qN, color = "darkblue") +
  xlab("") + ylab("") + ggtitle("loglikelihood") +theme_minimal() 

multiplot(h_plot, likelihood_plot,cols=2)





############j5###########################################################
x5 = read.table('X_j5.txt')
y5 = read.table('Y_j5.txt')
H_5 = read.table('output_H_vN_j5.txt')
#H_5 = read.table('output_H_qN_j5.txt')
#H_5 = read.table("output_H_gd_j5.txt")
data5 = data.frame("int" = x5[,1], "x" = x5[,2], "y" = y5[,1], "id" = y5[,2], "h1" = H_5[,1], "h2" = H_5[,2],"h3" = H_5[,3], "h4" = H_5[,4],"h5" = H_5[,5])
ggplot(data5, aes(x = x, y = y )) + geom_point()


maxiter = 2000
est_pars_5 = read.table('output_log_vN_j5.txt')
est_pars_5 = est_pars_5[1:maxiter,]

est_pars_5_qN = read.table('output_log_qN_j5.txt')
est_pars_5_qN = est_pars_5_qNR[1:maxiter,]


est_pars_5_gd = read.table('output_log_gd_j5.txt', skip = 50)
est_pars_5_gd = est_pars_5_gd[1:maxiter,]


idx_1_gd = 4
idx_2_gd = 3
idx_3_gd = 1
idx_4_gd = 2
idx_5_gd = 5

idx_1_qNR = 5
idx_2_qNR = 3
idx_3_qNR = 1
idx_4_qNR = 2
idx_5_qNR = 4


idx_1 = 5
idx_2 = 3
idx_3 = 4
idx_4 = 2
idx_5 = 1
p1 = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_1]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_1_gd]),  alpha = 0.1,data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_1_qNR]), alpha = 0.1, color = "darkblue") +
  geom_hline(yintercept = 2.5, color="red") + 
  xlab("iter") + ylab("expert1") + ggtitle("beta0") + theme(legend.position="none")

p2 = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_2]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_2_gd], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_2_qNR]), alpha = 0.1, color = "darkblue") +
  geom_hline(yintercept = 0.4, color="red") + 
  xlab("iter") + ylab("expert2")  + theme(legend.position="none")

p3 = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_3]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_3_gd], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_3_qNR]), alpha = 0.1, color = "darkblue") +
  geom_hline(yintercept = 1.2, color="red") + 
  xlab("iter") + ylab("expert3")  + theme(legend.position="none")

p4 = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_4]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_4_gd], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_4_qNR]), alpha = 0.1, color = "darkblue") +
  geom_hline(yintercept = -1.2, color="red") + 
  xlab("iter") + ylab("expert4")  + theme(legend.position="none")


p5 = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_5]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_5_gd], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_5_qNR]), alpha = 0.1, color = "darkblue") +
  geom_hline(yintercept = -1.4, color="red") + 
  xlab("iter") + ylab("expert5")  + theme(legend.position="none")

p6  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_1 + 5]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_1_gd + 5], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_1_qNR + 5], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  geom_hline(yintercept = -0.5, color="red") + 
  xlab("iter") + ylab("") + ggtitle("beta1") + theme(legend.position="none")

p7  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_2 + 5]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_2_gd + 5], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_2_qNR + 5], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  geom_hline(yintercept = 0.8, color="red") + 
  xlab("iter") + ylab("") +  theme(legend.position="none")

p8  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_3 +5]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_3_gd + 5], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_3_qNR + 5], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  geom_hline(yintercept = -0.1, color="red") + 
  xlab("iter") + ylab("") +  theme(legend.position="none")

p9  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_4 + 5]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_4_gd +5], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_4_qNR + 5], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  geom_hline(yintercept = 0.1, color="red") + 
  xlab("iter") + ylab("") +  theme(legend.position="none")

p10  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_5 +5]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_5_gd +5], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_5_qNR + 5], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  geom_hline(yintercept = 0.4, color="red") + 
  xlab("iter") + ylab("") +  theme(legend.position="none")

p11  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_1 +10]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_1_gd +10], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_1_qNR + 10], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  xlab("iter") + ylab("") + ggtitle("theta0")  +theme(legend.position="none")

p12  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_2 +10]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_2_gd +10], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_2_qNR + 10], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  xlab("iter") + ylab("")  +theme(legend.position="none")


p13  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_3 +10]),alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_3_gd +10], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_3_qNR + 10], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  xlab("iter") + ylab("")  +theme(legend.position="none")


p14  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_4 +10]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_4_gd +10], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_4_qNR + 10], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  xlab("iter") + ylab("")  +theme(legend.position="none")


p15  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_5 +10]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_5_gd +10], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_5_qNR + 10], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  xlab("iter") + ylab("")  +theme(legend.position="none")

p16  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_1 +15]),alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_1_gd +15], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_1_qNR + 15], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  xlab("iter") + ylab("") + ggtitle("theta1")  +theme(legend.position="none")

p17  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_2 +15]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_2_gd +15], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_2_qNR + 15], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  xlab("iter") + ylab("")   +theme(legend.position="none")

p18  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_3 +15]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_3_gd +15], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_3_qNR + 15], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  xlab("iter") + ylab("")  +theme(legend.position="none")

p19  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_4 +15]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_4_gd +15], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_4_qNR + 15], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  xlab("iter") + ylab("")  +theme(legend.position="none")

p20  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_5 +15]),alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_5_gd +15], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_5_qNR + 15], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  xlab("iter") + ylab("")  +theme(legend.position="none")

p21  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_1 + 20]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_1_gd +20], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_1_qNR + 20], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  geom_hline(yintercept = 0.09, color="red") + 
  xlab("iter") + ylab("") + ggtitle("sigma2")  +theme(legend.position="none")

p22  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_2 + 20]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_2_gd +20], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_2_qNR + 20], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  geom_hline(yintercept = 0.09, color="red") + 
  xlab("iter") + ylab("")  +theme(legend.position="none")

p23  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_3 + 20]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_3_gd +20], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_3_qNR + 20], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  geom_hline(yintercept = 0.09, color="red") + 
  xlab("iter") + ylab("")  +theme(legend.position="none")

p24  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_4 + 20]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_4_gd +20], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_4_qNR + 20], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  geom_hline(yintercept = 0.09, color="red") + 
  xlab("iter") + ylab("")  +theme(legend.position="none")

p25  = ggplot(data = est_pars_5) + theme_minimal() +  
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,idx_5 + 20]), alpha = 0.1, color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,idx_5_gd +20], alpha = 0.1), data = est_pars_5_gd, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qNR[,idx_5_qNR + 20], alpha = 0.1), data = est_pars_5_gd,  color = "darkblue") +
  geom_hline(yintercept = 0.09, color="red") + 
  xlab("iter") + ylab("")  +theme(legend.position="none")

multiplot(p1, p2, p3, p4, p5,p6,p7,p8,p9,p10, p11, p12, p13, p14, p15,p16,p17,p18,p19,p20,p21, p22, p23, p24, p25,cols=5)

h5_plot = ggplot(data5, aes(x, y = value,color = variable)) + 
  geom_point(aes(y = h1, col = "h1")) + 
  geom_point(aes(y = h2, col = "h2"))+ 
  geom_point(aes(y = h3, col = "h3"))+ 
  geom_point(aes(y = h4, col = "h4"))+ 
  geom_point(aes(y = h5, col = "h5"))+
  ggtitle("H") +theme_minimal()

colnames(est_pars_5_gd)[26] = "Gradient Descent"
colnames(est_pars_5_qN)[26] = "quasi-Newton"
colnames(est_pars_5)[26] = "variant Newton's"
likelihood5_plot = ggplot(data = est_pars_5_gd)  +theme_minimal()+
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_gd[,26]), alpha = 0.1, color = "grey") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5[,26], alpha = 0.1), data = est_pars_5,  color = "darkred") +
  geom_point(aes(1:nrow(est_pars_5), est_pars_5_qN[,26]), data = est_pars_5_qN, alpha = 0.1, color = "darkblue") +
  xlab("") + ylab("") + ggtitle("loglikelihood") 

multiplot(h5_plot, likelihood5_plot,cols=2)



