# 2018 11 26 I.Zliobaite
# analyze domestication

input_all <- 'data_input/candidates_traits.csv'
out_file <- 'outputs/predictions.csv'
out_file2 <- 'outputs/precision.csv'

set.seed(1982) # fixing randomization seed such that cross-validation is reproducible

no_cols <- 17

#read data
data_all <- read.csv(input_all, header = TRUE, sep = '\t')
colnames(data_all)[1] <- 'Genus'

# plot correlation between input variables
library(corrplot)
M <- cor(data_all[,2:no_cols])
pdf('plots/fig_correlations.pdf',width = 6, height = 6)
corrplot(M, method="circle")
dev.off()

# data visualization (PCA)
myPCA <- prcomp(data_all[,3:no_cols], scale = T)
pdf('plots/fig_pca.pdf',height = 11, width = 11)
ind <- which(data_all[,2]==1)
plot(myPCA$x[,1],myPCA$x[,2],col='white',xaxt='n',yaxt='n',ann=FALSE)
text(myPCA$x[,1],myPCA$x[,2],data_all[,1],cex=1.2, col = 'black')
text(myPCA$x[ind,1],myPCA$x[ind,2],data_all[ind,1],cex=1.2,col = 'red')
dev.off()

# fitting decision trees
library(rpart)
fit_tree <- rpart(Domesticated ~ ., method="class", data=data_all[,2:no_cols], minsplit = 10)
fit_tree2 <- rpart(Domesticated ~ ., method="class", data=data_all[,c(2:5,11:no_cols)], minsplit = 10)

# plotting trees
pdf('plots/fig_tree.pdf')
plot(fit_tree, uniform=TRUE) 
text(fit_tree, use.n=TRUE, all=TRUE, cex=.8)
dev.off()

pdf('plots/fig_tree2.pdf')
plot(fit_tree2, uniform=TRUE) 
text(fit_tree2, use.n=TRUE, all=TRUE, cex=.8)
dev.off()

# probability estimates for domestication using trees
pred_tree <- predict(fit_tree, data_all[,2:no_cols],type = "prob")
pred_tree <- round(pred_tree,digits = 2)
pred_tree_cl <- predict(fit_tree, data_all[,2:no_cols],type = "class")

pred_tree2 <- predict(fit_tree2, data_all[,2:no_cols],type = "prob")
pred_tree2 <- round(pred_tree2,digits = 2)
pred_tree2_cl <- predict(fit_tree2, data_all[,2:no_cols],type = "class")


# fiting regression
library(glmnet)
x <- as.matrix(data_all[,3:no_cols])
y <- as.vector(data_all[,2])
cv.out <- cv.glmnet(x,y,alpha=0.5,family='binomial',type.measure = 'deviance',nfolds = 10)
pdf('plots/fig_cross_validation.pdf')
plot(cv.out)
dev.off()

#min value of lambda
lambda_min <- cv.out$lambda.min
#best value of lambda
lambda_1se <- cv.out$lambda.1se
#regression coefficients
print('resulting logistic regression (all variables) model coefficients:')
print(coef(cv.out,s=lambda_min))

# probability estimates for domestication using regression
pred_las <- predict(cv.out, x, s="lambda.min",type = 'response')
pred_las <- round(pred_las,digits = 2)


# regression only ecological
x2 <- as.matrix(data_all[,c(3:9,11:no_cols)])
cv.out2 <- cv.glmnet(x2,y,alpha=0.5,family='binomial',type.measure = 'deviance',nfolds = 10)
pdf('plots/fig_cross_validation2.pdf')
plot(cv.out2)
dev.off()

#min value of lambda
lambda_min2 <- cv.out2$lambda.min
#best value of lambda
lambda_1se2 <- cv.out2$lambda.1se
#regression coefficients
print('resulting logistic regression (no "tropical") model coefficients:')
print(coef(cv.out2,s=lambda_min2))

pred_las2 <- predict(cv.out2, x2, s="lambda.min",type = 'response')
pred_las2 <- round(pred_las2,digits = 2)


# recording predictions
data_out <- cbind(data_all[,c('Genus','Domesticated')],pred_las,pred_las2,pred_tree[,2],pred_tree2[,2])
colnames(data_out)[3:6] <- c('predictions_regression','predictions_regression2','predictions_tree_all','predictions_tree_org')
ord <- order(-data_out[,'predictions_regression2'])
data_out <- data_out[ord,]
write.table(data_out,file = out_file,row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE )

# compute mean suitability scores and topx rates

print('Standard deviations of suitability scores accross realms')

data_rates <- c()
realms <- c('Nearctic','Palearctic','Afrotropic','Indomalaya','Australasia','Neotropic')
nn <- dim(data_all)[1]

for (sk in 0:nn){
  if (sk<=1){
    topx <- rep(0,nn)
  }
  if (sk>0){
    thresholds <- pred_las2[order(-pred_las2)]
    ind_now <- which(pred_las2 == thresholds[sk])
    n_now <- length(ind_now)
    if (n_now==1){
      topx[ind_now] <- 1
    }else{
      sk_back <- length(which(pred_las2 >= thresholds[sk])) - n_now
      topx[ind_now] <- (sk - sk_back)/n_now
    }
    if (sum(topx)!=sk){
      print('problem 1')
    }
  }else{
    topx <- pred_las2
  }
  vec_realms <- c()
    for (rr in realms){
      ind <- which(!is.na(data_all[,rr]))
      rate_est <- 100*sum(topx[ind]*data_all[ind,rr])/ sum(data_all[ind,rr])
      if (sk==0){
        print(rr)
        print(sd(topx[ind]))
      }
      vec_realms <- c(vec_realms,rate_est)
    }
  data_rates <- rbind(data_rates,c(sk,round(vec_realms,digits = 2)))
}
colnames(data_rates) <- c('k',realms)
write.table(data_rates,file = out_file2,row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE )

pdf('plots/fig_precision.pdf', width = 6.5, height = 5)
rg <- 2:69
lw <- 4
ord_cont <- c(2,1,6,3,4,5)
mycol <- c("#7fc97f","#000000", "#f0027f", "#386cb0", "#ffff99", "#fdc086")
plot(NA,NA,xlim = c(1,68),ylim = c(0,100),ylab = 'Potential domestication rate, %', xlab = 'Number of best candidates to domesticate')
lines(c(10,10),c(0,100),col = 'grey',type = 'l')
lines(data_rates[rg,1],data_rates[rg,7],type='l',col = mycol[6], lwd = lw)
lines(data_rates[rg,1],data_rates[rg,6],type='l',col = mycol[5], lwd = lw)
lines(data_rates[rg,1],data_rates[rg,5],type='l',col = mycol[4], lwd = lw)
lines(data_rates[rg,1],data_rates[rg,4],type='l',col = mycol[3], lwd = lw)
lines(data_rates[rg,1],data_rates[rg,2],type='l',col = mycol[1], lwd = lw)
lines(data_rates[rg,1],data_rates[rg,3],type='l',col = mycol[2], lwd = lw)
points(10,27.27,col = mycol[2], cex = 1.5,lwd = 2)
points(10,14.28,col = mycol[4], cex = 1.5,lwd = 2)
points(10,28.57,col = mycol[6], cex = 1.5,lwd = 2)
points(10.3,0,col = mycol[5], cex = 1.5,lwd = 2)
points(9.7,0,col = mycol[3], cex = 1.5,lwd = 2)
points(10,0,col = mycol[1], cex = 1.5,lwd = 2)
legend(12,100,realms[ord_cont],lty = 1,col = mycol[ord_cont], cex = 0.75, lwd = lw, bty = "n")
dev.off()

