
#R CMD INSTALL -l /home/kjlee/R/x86_64-redhat-linux-gnu-library/3.3 /data8/kjlee/Research/MixtureRegressionModels/Kuo-Jung/Code/Package/UnknownCompFMR_2.0.tar.gz

require(MASS)
require(pscl)
require(msm)
#require(flexmix)
require(gtools)
rm(list=ls(all=TRUE))

NORMALIZATION = TRUE
Robust = 1
RJMCMC = 1


load("/data8/kjlee/Research/MixtureRegressionModels/Yi-Chi/Data/CreditGrowth_JIMF.rda")


#makes dataset balanced, only complete data are analysed
data.num = 1
data.name = colnames(X)[data.num]

y = (X[, data.num])

#y = sqrt(y+2)
x = as.matrix(X[, -data.num])

cov.of.no.interest = NULL
cov.of.interest = colnames(x)[!colnames(x)%in%cov.of.no.interest]

x = x[, colnames(x)%in%cov.of.interest]


#Determine whether normalizing the design matrix, X. Dummy varibles are not be normalized.
if(NORMALIZATION){
#y = (y-mean(y))

dummy.variable = c("oil.exporter", "floater", "infl.targeter")

non.dummy.x = x[, !(colnames(x) %in% dummy.variable)]


x.normalized = apply(non.dummy.x, 2, function(xx) (xx-mean(xx))/sd(xx))

x[, !(colnames(x) %in% dummy.variable)] = x.normalized
}

x = cbind(1, x)
colnames(x) = c("Intercept", cov.of.interest)
#Whether Intercept is included. Intercept is always included as the paper suggested.


#==========================================================================#
num.of.iteration= 10000
num.of.iterations.Inside = 10
#==========================================================================#
max.num.of.groups = 5
group.info = rep(NA, max.num.of.groups) 

#==========================================================================#

#num.of.groups = G = 2

#crisis = cbind(y, x)
#crisis = as.data.frame(crisis)
#flexmix.fit = flexmix(y~., data=crisis, k=2, control=list(verb=0, iter=100))
#beta.est = cbind(parameters(flexmix.fit, component = 1), parameters(flexmix.fit, component = 2))
#beta.est = t(beta.est[-c(1, nrow(beta.est)), ])

initial.num.of.group = 1

num.of.covariates =  ncol(x)
num.of.obs = length(y)
num.of.obs.for.each.group = round(num.of.obs*rep(1/initial.num.of.group, initial.num.of.group))	
num.of.obs.for.each.group[initial.num.of.group] = num.of.obs - sum(num.of.obs.for.each.group[1:(initial.num.of.group-1)])


c2 = 100
alpha.g = rep(1, max.num.of.groups)

sigma2.a = sigma2.b  = 1
dG = bG = 0.5
nu = 4
tuning.para = 0.05
dgj = rep(0.5, num.of.covariates) #c(rep(0.9, 12), rep(0.01, num.of.covariates-12))

z =  sample(0:(initial.num.of.group-1), num.of.obs, replace=TRUE)

#==================================================================================#

r.samples = array(NA, c(max.num.of.groups, (num.of.covariates), num.of.iteration))
r.sample = matrix(r.samples[, , 1], nrow=max.num.of.groups)
r.sample[1:initial.num.of.group, ] = rbinom(length(r.sample[1:initial.num.of.group, ]), 1, 0.5)
r.samples[1:initial.num.of.group,, 1] = r.sample[1:initial.num.of.group, ]

rho.samples = matrix(NA, max.num.of.groups, num.of.iteration)
rho.sample = rho.samples[, 1]
rho.sample[1:initial.num.of.group]= 1/initial.num.of.group
rho.samples[, 1] = rho.sample 

omega.samples = matrix(NA, num.of.obs, num.of.iteration)
omega.sample = rep(1, num.of.obs) 
omega.samples[, 1] = omega.sample
#group.info[as.numeric(dimnames(table(z))[[1]])+1] = as.vector(table(z))

#==================================================================================#

InitialValues = list(z = z, r.sample = r.sample, rho.sample = rho.sample, omega.sample = omega.sample)
Data = list(Y=y, X=x)
GivenValues = list(alpha.g = alpha.g, c2 = c2, sigma2.a = sigma2.a, sigma2.b = sigma2.b, nu=nu, 
	tuning.para = tuning.para, dgj = dgj)

#install.packages("UnknownCompFMR_1.0.tar.gz")
library(UnknownCompFMR)
posterior.samples = VariableSelectionUnknowCompFMR(num.of.iteration, num.of.iterations.Inside, max.num.of.groups, Data, 
	InitialValues, GivenValues, Robust, RJMCMC)

save(posterior.samples, file="CreditGrowth_PosteriorSamples_Demean_InitialGroup5_Inner_10_c2_100.RData")


print(table(apply(posterior.samples$z, 2, function(x) length(unique(x)))))

#==================================================================================#

if(0){
num.of.groups = 1
z.samples = posterior.samples$z
r.samples = posterior.samples$r.samples
rho.samples = posterior.samples$rho.samples

group.num.list = apply(z.samples, 2, function(x) length(unique(x)))

z.samples.LS = z.samples[, group.num.list == num.of.groups] 
r.samples.LS = r.samples[,,group.num.list== num.of.groups]
rho.samples.LS = rho.samples[, group.num.list == num.of.groups]

num.of.iter.eff = ncol(z.samples.LS)

r.samples.rm.nan = array(0, c(num.of.groups, num.of.covariates, num.of.iter.eff))


num.of.obs.in.groups = apply(z.samples.LS, 2, table)

m.keep = 0
for(m in 1:num.of.iter.eff){
	r.tmp = r.samples.LS[,,m][complete.cases(r.samples.LS[,,m]*0), ,drop=FALSE]
	if(nrow(r.tmp) == num.of.groups ){
		m.keep = c(m.keep, m)
		r.samples.rm.nan[, , m] = r.samples.LS[,,m][complete.cases(r.samples.LS[,,m]*0), ,drop=FALSE]
		if(num.of.obs.in.groups[1]<num.of.obs.in.groups[2])
		r.samples.rm.nan[, , m] = r.samples.rm.nan[c(2,1), , m] 
	}
}


r.est = apply(r.samples.rm.nan[,,m.keep, drop=FALSE], c(1, 2), mean)
colnames(r.est) = colnames(x)

r.est = as.table(r.est)

pdf(file="gammabarchart2.pdf", height = 14, width = 7)
par(las=2)
par(mar=c(5,14,4,2))

barplot(r.est, main= expression(paste("Posterior estimate of P(", gamma, "=1)")), 
xlab="Probability", col=c("red", "blue"),
	names.arg = colnames(r.est), beside=TRUE, horiz=TRUE)


#barplot(r.est, main= expression(paste("Posterior estimate of P(", gamma, "=1)")), 
#  xlab="Probability", col=c("red"), names.arg = colnames(r.est), beside=TRUE, horiz=TRUE)

abline(v=0.5, lwd=2, col=2)

dev.off()
}

#print(sum(posterior.samples$z[, 1] == posterior.samples$z[, num.of.iteration])/num.of.obs)
#print(apply(posterior.samples$r.samples, c(1, 2), mean))
#print(apply(posterior.samples$rho.samples, 1, mean))

