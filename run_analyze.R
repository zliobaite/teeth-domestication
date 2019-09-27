# 2018 11 26 I.Zliobaite
# analyze domestication

input_all <- 'data_input/candidates_traits.csv'
out_file <- 'outputs/predictions.csv'

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
print('resulting logistic regression model coefficients:')
print(coef(cv.out,s=lambda_min))

# probability estimates for domestication using regression
pred_las <- predict(cv.out, x, s="lambda.min",type = 'response')
pred_las <- round(pred_las,digits = 2)

# compute mean suitability scores
realms <- c('Nearctic','Palearctic','Afrotropic','Indomalaya','Australasia','Neotropic')
print('Mean suitability scores accross realms')
for (rr in realms){
  ind <- which(!is.na(data_all[,rr]))
  print(rr)
  print('regression')
  print(round(sum(pred_las[ind]*data_all[ind,rr])/ sum(data_all[ind,rr]),digits = 2))
}


# recording predictions
data_out <- cbind(data_all[,c('Genus','Domesticated')],pred_las,pred_tree[,2],pred_tree2[,2])
colnames(data_out)[3:5] <- c('predictions_regression','predictions_tree_all','predictions_tree_org')
ord <- order(-data_out[,'predictions_regression'])
data_out <- data_out[ord,]
write.table(data_out,file = out_file,row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE )

