# Check the speed of semi-parallel approach
# for logistic regression with covariates
# Karolina Sikorska and Paul Eilers, 2013

#####sample dataset from the paper#########
#set.seed(2013)
#n = 10000
#m = 1000
#k = 1
#S = matrix(2 * runif(n * m), n, m)
#X0 = matrix(rnorm(n * k), n, k)
#X = cbind(1, X0)
#y = rbinom(n, size = 1, prob = c(0.5, 0.5))
#####sample dataset from the paper##########


######import dataset extracted from PGP########
#SNPs
S_list = read.table("snpMat.txt", header=T)
S = matrix(unlist(S_list, use.names=F), ncol = length(S_list), byrow=F)

covariate = read.csv("covariates.csv", header=T, sep=",")

#disease/control condition
y = unlist(covariate[2], use.names=F)

#covariates w/ missing values
rawX0 = covariate[3:5]

#norm
normalize <- function(x){(x-min(x))/(max(x)-min(x))}


#impute missing values w/ mean
X0_list = lapply(rawX0, function(x) ifelse(is.na(x), mean(x, na.rm=T), x))
X0 = (matrix(unlist(X0_list), ncol = length(X0_list), byrow=F))
X0 = normalize(X0)




#covariates
X = cbind(1, X0)

#number of individuals
n = nrow(S)
#number of SNPs
m = ncol(S)
#number of covariates
k = ncol(X0)
######import dataset extracted from PGP########

cat(sprintf("the number of samples : %d \n", n))
cat(sprintf("the number of SNPs : %d \n", m))
cat(sprintf("the number of covariates : %d \n", k))

# Do the computations
t0 = proc.time()[1]
mod0 = glm( y ~ X-1, family = binomial("logit")) 
p = mod0$fitted
w = p * (1 - p)
z = log(p / (1 - p)) + (y - p) / (p * (1 - p))
xtw = t(X * w)
U1 = xtw %*% z
U2 = solve(xtw %*% X, U1)
ztr = z  - X %*% U2
U3 = xtw %*% S
U4 = solve(xtw %*% X, U3)
Str = S - X %*% U4
Str2 = colSums(w * Str^2)
b = crossprod(ztr * w, Str)/Str2
err = sqrt(1/ Str2)
pvalue_parallel = 2 * pnorm(-abs(b / err))
pval = as.numeric(as.matrix(pvalue_parallel))

statistics_parallel = (b * b) / (err * err)
tstat = t(statistics_parallel)
tpval = t(pvalue_parallel)
bind = cbind(tpval, tstat)
bind = cbind(bind, errorsq)
#statistics = as.numeric(as.matrix(statistics_parallel))
errorsq = err * err

# Report time
t1 = proc.time()[1] - t0
msip = 1e-06 * n * m / t1
cat(sprintf("Speed: %2.1f Msips\n", msip))

###check results

t2 = proc.time()[1]
b1 = NULL
err1 = NULL
pval1 = NULL
for(i in 1:m){
  mod1 = glm(y ~ S[, i] + X0, family = binomial ("logit"))  
  b1[i] = summary(mod1)$coef[2, 1]
  err1[i] = summary(mod1)$coef[2, 2]
  pval1[i] = summary(mod1)$coef[2, 4]
}


t3 = proc.time()[1] - t2
msip1 = 1e-06 * n * m / t3
cat(sprintf("Speed: %2.1f Msips\n", msip1))



#true positive
print("true positive rate")
length(which(pval1[which(pval<=0.01)]<=0.01))/length(which(pval<=0.01))


#type i error/false positive
print("type i error/false positive rate")
length(which(pval1[which(pval>0.01)]<=0.01))/(length(which(pval1[which(pval>0.01)]<=0.01)) + length(which(pval1[which(pval>0.01)]>0.01)))


#type ii error/false negative
print("type ii error/false negative rate")
length(which(pval1[which(pval<=0.01)]>0.01))/length(which(pval<=0.01))

#true negative
print("true negative rate")
length(which(pval1[which(pval>0.01)]>0.01))/length(which(pval>0.01))

#test
c = 7
xtw1 = X * 0.25
inv1 = solve(t(xtw1) %*% X, t(X))
ddd1 = det(t(xtw1) %*% X)
beta1 = c * ddd1 * inv1 %*% (y - 0.5)
p1 = 1 / (1 + exp(-X%*%beta1))
w1 = p1 * (1 - p1)
test1 = X %*% beta1

xtw2 = X * w1[,1]
inv2 = solve(t(xtw2) %*% X, t(X))
ddd2 = det(t(xtw2) %*% X)
beta2 = beta1 + c * ddd2 * inv2 %*% (y - p1)
p2 = 1 / (1 + exp(-X%*%beta2))
w2 = p2 * (1 - p2)
test2 = X %*% beta2

xtw3 = X * w2[,1]
inv3 = solve(t(xtw3) %*% X, t(X))
ddd3 = det(t(xtw3) %*% X)
beta3 = beta2 + c * ddd3 * inv3 %*% (y - p2)
p3 = 1 / (1 + exp(-X%*%beta3))
w3 = p3 * (1 - p3)
test3 = X %*% beta3

xtw4 = X * w3[,1]
inv4 = solve(t(xtw4) %*% X, t(X))
ddd4 = det(t(xtw4) %*% X)
beta4 = beta3 + c * ddd4 * inv4 %*% (y - p3)
p4 = 1 / (1 + exp(-X%*%beta4))
w4 = p4 * (1 - p4)
test4 = X %*% beta4

xtw5 = X * w4[,1]
inv5 = solve(t(xtw5) %*% X, t(X))
ddd5 = det(t(xtw5) %*% X)
beta5 = beta4 + c * ddd5 * inv5 %*% (y - p4)
p5 = 1 / (1 + exp(-X%*%beta5))
w5 = p5 * (1 - p5)
test5 = X %*% beta5

xtw6 = X * w5[,1]
inv6 = solve(t(xtw6) %*% X, t(X))
ddd6 = det(t(xtw6) %*% X)
beta6 = beta5 + c * ddd6 * inv6 %*% (y - p5)
p6 = 1 / (1 + exp(-X%*%beta6))
w6 = p6 * (1 - p6)
test6 = X %*% beta6

xtw7 = X * w6[,1]
inv7 = solve(t(xtw7) %*% X, t(X))
ddd7 = det(t(xtw7) %*% X)
beta7 = beta6 + c * ddd7 * inv7 %*% (y - p6)
p7 = 1 / (1 + exp(-X%*%beta7))
w7 = p7 * (1 - p7)
test7 = X %*% beta7

xtw8 = X * w7[,1]
inv8 = solve(t(xtw8) %*% X, t(X))
ddd8 = det(t(xtw8) %*% X)
beta8 = beta7 + c * ddd8 * inv8 %*% (y - p7)
p8 = 1 / (1 + exp(-X%*%beta8))
w8 = p8 * (1 - p8)
test8 = X %*% beta8