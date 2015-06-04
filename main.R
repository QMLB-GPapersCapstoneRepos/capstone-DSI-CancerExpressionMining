setwd("C:/Users/Shin/Google Drive/SYS 6016/final_project")

source("http://bioconductor.org/biocLite.R")
# http://www.bioconductor.org/packages/release/bioc/html/made4.html
#biocLite("made4")
#biocLite("made4")

library("openxlsx")
library("made4")
library("ade4")
library(glmnet)
library(mht)
library(CePa)
library(made4)
library(ade4)

source("bolasso.R")


# read data, name of gene as row, name of patient as column
mydf <- read.xlsx("5998+array.xlsx", sheet = 1, startRow = 3, colNames = TRUE)

# Transform data by transposing so that name of gene as column, and name of patient
# as row. 45 observations in total. (45 patient: 8 normal people, 37 infected patients)
Tmydf<-t(mydf[,c(-1,-2)])
colnames(Tmydf) <- mydf[,1]
Tmydf<-data.frame(Tmydf)



# Check the missing value
missing <- is.na(Tmydf)
which(missing == TRUE)
# integer(0)
# No missing value exist

####### Create heatmap with respect to the original data values #######
# calculate Variance, ncol(Tmydf) = 10443
result<-c()               
for (i in (1:ncol(Tmydf))){
  result[i]<-var(Tmydf[,i])
}

# sort variance calculated per gene and get top 100 largest variance
Sortresult<-sort(result,decreasing=T)

# corresponding postion 
Top100postion<-match(Sortresult[1:100],result)

#How much total variance the top 100 genes capture
sum(Sortresult[1:100])/sum(Sortresult)
#about 41.91025% of total variance is captured by top 100 genes.

# get the dataset with genes having the top 100 largest variance
Top100dataset<-Tmydf[,c(Top100postion)]
m_matrix<-data.matrix(Top100dataset)

# heatmap #1.2
## get pearson correlation and clustering 
cor_t <- cor(t(m_matrix))
distancet <- as.dist((1-cor_t)/2)
hclust_complete <- hclust(distancet, method = "complete")
dendcomplete <- as.dendrogram(hclust_complete)

# plot heatmap 
heatmap(m_matrix, Rowv=dendcomplete, Colv=NA, scale="column")

## heatmap with 2-way hierarchical clustering
heatmap.2(t(m_matrix), trace="none", density="none", 
          scale="row",
          labRow="",
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2))

###### Create heatmap with respect to standardize for top 100 genes #######

# Variance calculation, ncol(Tmydf.sdize) = 10443              
result<-c()               
for (i in (1:ncol(Tmydf))){
  result[i]<-var(Tmydf[,i])
}

Sortresult<-sort(result,decreasing=T)

Top100postion<-match(Sortresult[1:100],result)
Top100dataset<-Tmydf[,c(Top100postion)]

# standardize the dataset
standardize<-function(x){return((x-mean(x))/(sd(x)))}
Top100dataset.sdize<-as.data.frame(lapply(Top100dataset,standardize))
rownames(Top100dataset.sdize)<-rownames(Tmydf)
summary(Top100dataset.sdize)


# heatmap #2.2
## pearson
cor_t <- cor(t(m_matrix.norm))
distancet <- as.dist((1-cor_t)/2)
hclust_complete <- hclust(distancet, method = "complete")
dendcomplete <- as.dendrogram(hclust_complete)
heatmap(as.matrix(Top100dataset.sdize), Rowv=dendcomplete, Colv=NA, scale="column")

## 2-way
heatmap.2(t(as.matrix(Top100dataset.sdize)), trace="none", density="none", 
          scale="row",
          labRow="",
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2))
####### End heatmap #######



######################## lasso regression ####################################
# standardize
standardize<-function(x){return((x-mean(x))/(sd(x)))}
Tmydf.sdize<-as.data.frame(lapply(Tmydf,standardize))


# Add an attribute "class" to distinguish normal people and patients, normal people as 0 and patients as 1.
Tmydf.sdize$class<-1
#Frist 8 normal people as 0;
Tmydf.sdize$class[c(1:8)]<-0
summary(Tmydf.sdize)


#we utilized standardized data set here
y <- as.numeric(Tmydf.sdize[,10444])
y.f <- as.factor(y)
x <- data.matrix(Tmydf.sdize[,1:10443])
fit = glmnet(x, y, family = "binomial",standardize = F)
## Coefficient path plot
plot(fit,xvar = "lambda",label = TRUE)

##Choose optimal lambda via Misclassification Error
cvfit = cv.glmnet(x, y, nfolds=nrow(x),family = "binomial", type.measure = "class",standardize = F)
# Misclassification Error plot
plot(cvfit)


#value of lambda that gives minimum AUC.
cvfit$lambda.min
#largest value of lambda such that error is within 1 standard error of the minimum.
#Normally, people use lambda.1se for model simplex 
cvfit$lambda.1se

# Coeffecient Plot
plot(fit,xvar = "lambda",label = TRUE)
abline(v=log(cvfit$lambda.min), lty=2)
abline(v=log(cvfit$lambda.1se), lty=2)
abline(h=0, lty=2)
coef<-coef(cvfit, s = "lambda.1se")
str(coef)
summary(coef)

#get coef not equal to 0
a<-c()
for (i in 1:(ncol(Tmydf)-1)){
  a[i]<-coef[i,]!=0
}

#get the position 
b<-which(a=="TRUE")

#since first one is interception 
c=b-1

#get the gene name
d<-Tmydf.sdize[,c]
colnames(d)

# get the value of coefficients
coef.result <- data.frame(colnames(d),coef@x[2:length(b)])

# sort result by absolute value
coef.result.sorted <- coef.result[order(abs(coef.result[,2]), decreasing = T),]
View(coef.result.sorted)


# bootstrap lasso ----


# Order variables using bolasso
dyadiqueordre(x, y, m = 100, maxordre = 45, family = "binomial")


# SPINK2  GOLGA8A  LYZ  NAGPA  DUSP2  NLGN2  CD34  CASK  PHGDH  NELL2  FMNL3  HIPK1  SERF1A  SPRY1  CBWD3  
# PAPD1  ACSM3  MAP6D1  TXK  LAPTM4B  NR4A3  MIR573  CAMK2N1  SNORA4  TNFSF9  MIA3  AAK1  NEK3  GNRH1  ZNF331 
# ZSWIM4 INADL DOCK3 IFT74 LMBRD2 SNORA31 SNORA27 MGC18216 HPR NAPSB SORL1 DHX34 LOC644936 GZMK SYNJ2BP B4GALT7 BIN3

#####################################################
############ rule fit ###############################
####################################################
platform = "windows"
rfhome = "C:/Users/Shin/Documents/R/win-library/3.1/RuleFit3"
source("C:/Users/Shin/Documents/R/win-library/3.1/RuleFit3/rulefit.r")
#install.packages("akima", lib=rfhome)
library(akima, lib.loc=rfhome)

y.n <- as.numeric(y)
y.n[y.n==0] <- -1
x.n <- as.matrix(x) 
str(x.n)

rfmod <- rulefit(x.n,y.n, sparse = 1, rfmode="class", max.rules = 2000, tree.size =2)
summary(rfmod)
str(rfmod)


#par(mfrow=c(1,1))
vi = varimp(1:20, x = x.n)
runstats(rfmod)
rfmodinfo(rfmod)

modxval <- rfxval (nfold=5, quiet=F)
vi$imp[1:30]
vi = varimp(1:20, x = x.n)

bost.null= intnull ()
interact (vi$ord[1:20], bost.null)
interact (vi$ord[1:20])

int2var= twovarint("KIT", c("SPINK2", "LYZ", "SNORA4"), null.mods = bost.null, ymax = 0.02)
################### End of first data set ##################

#################### second data set ####################
mydy2<-read.gct("Capstone2.gct")
data.2<-mydy2
data.2<-data.frame(t(data.2)
                   
####set normal donor as 0 and patients as 1
data.2$class<-1
data.2$class[c(1:6)]<-0

# Check the missing value
missing <- is.na(data.2)
which(missing == TRUE)
# integer(0)
# No missing value exist


###################### Visualization ####################




# calculate Variance, ncol(Tmydf) = 10443
result<-c()               
for (i in (1:ncol(data.2))){
  result[i]<-var(data.2[,i])
}

# sort variance calculated per gene and get top 100 largest variance
Sortresult.2<-sort(result,decreasing=T)

# corresponding postion 
Top100postion.2<-match(Sortresult.2[1:100],result)

# get the dataset with genes having the top 100 largest variance
Top100dataset.2<-data.2[,c(Top100postion.2)]

m_matrix.2 <- data.matrix(Top100dataset.2)

summary(m_matrix.2)
# heatmap #2.2
## pearson correlation
cor_t <- cor(t(m_matrix.2))
distancet <- as.dist((1-cor_t)/2)
hclust_complete <- hclust(distancet, method = "complete")
dendcomplete <- as.dendrogram(hclust_complete)
heatmap(m_matrix.2, Rowv=dendcomplete, Colv=NA, scale="column")


## 2-way
heatmap.2(t(m_matrix.2), trace="none", density="none", 
          scale="row",
          labRow="",
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2))



# standardize across genes 
data.2.sdize<-as.data.frame(lapply(data.2,standardize))
rownames(data.2.sdize) <- rownames(data.2)


#Sdize.genes<-rank(data.2.sdize,100)
Top100.std.2<-data.2.sdize[,c(Top100postion.2)]

m_matrix.std <- as.matrix(Top100.std.2)
heatmap(m_matrix.std, Colv=NA, scale="column")


# # standardization across genes 
# standardize<-function(x){return((x-mean(x))/(sd(x)))}
# data.2.sdize<-as.data.frame(lapply(data.2,standardize))
# 
# 
# 
# Stan.gene<-rank(data.2.sdize,100)
# m_matrix1 <- data.matrix(Stan.gene)
# m_matrix1 <- data.matrix(t(Stan.gene))
# heatmap(m_matrix1, Colv=NA, scale="column")

# looks so wired!!!!!

###################### Lasso Regression ##########################

##standardize###
mydy2.sdize<-as.data.frame(lapply(mydy2,standardize))
summary(mydy2.sdize)

######
y2 <- as.numeric(mydy2.sdize[,20726])
x2 <- data.matrix(mydy2.sdize[,1:20725])

fit2 = glmnet(x2, y2)


plot(fit2,xvar = "lambda",label = TRUE)

#Choose optimal lambda
cvfit2<-cv.glmnet(x2,y2,nfolds=45,family = "binomial", type.measure = "class")
plot(cvfit2)


#value of lambda that gives minimum cvm.
cvfit2$lambda.min
#largest value of lambda such that error is within 1 standard error of the minimum.
#Normally, people use lambda.1se for model simplex 
cvfit2$lambda.1se

#final plot
plot(fit2,xvar = "lambda",label = TRUE)
abline(v=log(cvfit2$lambda.1se), lty=2)
abline(h=0, lty=2)

coef2<-coef(cvfit2, s = "lambda.1se")

#get coef not equal to 0
e<-c()
for (i in 1:10444){
  e[i]<-coef2[i,]!=0
}

#get the position 
f<-which(e=="TRUE")

#since first one is interception 
g=f-1

#get the gene name
h<-mydy2[,g]
colnames(h)

# bootstrap lasso -----

# Order variables using bolasso
dyadiqueordre(x.2, y.2, m = 100, maxordre = 45, family = "binomial")

