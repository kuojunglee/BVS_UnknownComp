
#R CMD INSTALL -l /home/kjlee/R/x86_64-redhat-linux-gnu-library/3.3 /data8/kjlee/Research/MixtureRegressionModels/Kuo-Jung/Code/Package/UnknownCompFMR_2.0.tar.gz

require(MASS)
require(pscl)
require(msm)
library(rgdal)
library(ggplot2)
library(maptools)
require(gdata)
require(gtools)
library(rworldmap)
library(gpclib)
library(flexmix)
library(UnknownCompFMR)

rm(list=ls(all=TRUE))
set.seed(1)

NORMALIZATION = TRUE
Robust = 1
RJMCMC = 1
ClusteringMap = 0

if(Sys.info()[1]!="Darwin"){
load("/home/kuojung/Research/KJLEE_Papers/MixtureRegressionModels/Yi-Chi/Data/CreditGrowth_JIMF.rda")
}
if(Sys.info()[1]=="Darwin"){
load("/Users/kjlee/Research/KJLEE_Papers/MixtureRegressionModels/Yi-Chi/Data/CreditGrowth_JIMF.rda")
}


#makes dataset balanced, only complete data are analysed
data.num = 1
data.name = colnames(X)[data.num]

y = (X[, data.num])

if(0){
CumulativeLoss.df = data.frame(CumulativeLoss = y)

ggplot(CumulativeLoss.df, aes(x=CumulativeLoss)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2) + 
 labs(x = "Cumulative Loss")+
ggsave("CULOSSBlack.png", plot = last_plot(), device = "png", path = "/Users/kjlee/Research/KJLEE_Papers/MixtureRegressionModels/YC_Martin_KJ/Figures/",
       scale = 0.55, width = 12, height = 10, units = c("in"), dpi = 300, limitsize = TRUE)
}
#y = sqrt(y+2)
x = as.matrix(X[, -data.num])

cov.of.no.interest = c("openness.0206", "cpi.corruption.06") #, "merchTrade.0006", "freedom.from.corr.06")
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
num.of.iteration= 1000
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

initial.num.of.group =5

num.of.covariates =  ncol(x)
num.of.obs = length(y)
num.of.obs.for.each.group = round(num.of.obs*rep(1/initial.num.of.group, initial.num.of.group))	
num.of.obs.for.each.group[initial.num.of.group] = num.of.obs - sum(num.of.obs.for.each.group[1:(initial.num.of.group-1)])


c2 = num.of.obs
lambda = 3
alpha.g = rep(1, max.num.of.groups)

sigma2.a = 10
sigma2.b  = 10
dG = bG = 0.5
nu = 4
tuning.para = 4 #using in MH to generate omega N(omega, tuning para)
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
	tuning.para = tuning.para, dgj = dgj, lambda = lambda)

#install.packages("UnknownCompFMR_1.0.tar.gz")
posterior.samples = VariableSelectionUnknowCompFMR(num.of.iteration, num.of.iterations.Inside, max.num.of.groups, Data, 
	InitialValues, GivenValues, Robust, RJMCMC)


print(table(apply(posterior.samples$z, 2, function(x) length(unique(x)))))

if(0){

load("/Users/kjlee/Research/KJLEE_Papers/MixtureRegressionModels/Kuo-Jung/Code/R_Code/CreditGrowth_PosteriorSamples_Demean_InitialGroup1_Inner_10_c2_100.RData")
num.of.groups = 2
z.samples = posterior.samples$z
r.samples = posterior.samples$r.samples
rho.samples = posterior.samples$rho.samples

#group.num.list = apply(z.samples, 2, function(x) length(unique(x)))

group.num.list = apply(r.samples[, 1,], 2, function(x) sum(is.finite(x)))

z.samples.LS = z.samples[, group.num.list == num.of.groups] 
r.samples.LS = r.samples[,,group.num.list == num.of.groups, drop=FALSE]
rho.samples.LS = rho.samples[, group.num.list == num.of.groups, drop=FALSE]

omega.samples = posterior.samples$omega.samples[, group.num.list == num.of.groups, drop=FALSE]

num.of.iter.eff = ncol(z.samples.LS)

num.of.iter.eff.remove = 0

r.samples.rm.nan = array(0, c(num.of.groups, num.of.covariates, (num.of.iter.eff-num.of.iter.eff.remove)))
rho.samples.rm.nan = matrix(0, num.of.groups, (num.of.iter.eff-num.of.iter.eff.remove))


num.of.obs.in.groups = apply(z.samples.LS, 2, table)



#=======================================================================#
	prob.cube = prob.cube.ICL = prob.cube.cluster =array(0, c(num.of.iter.eff, num.of.obs, num.of.groups))
	prob.cube.logLike = 0
	beta.est.group = array(0, c(num.of.iter.eff, num.of.covariates, num.of.groups))
	sigma.est.group = matrix(0, num.of.iter.eff, num.of.groups)
	for(m in 1:num.of.iter.eff){
		group.index = unique(z.samples.LS[, m])
		for(g in 1:num.of.groups){
			beta.est = matrix(0, ncol(r.samples.LS[,,m]), 1) 
			y.g = y[z.samples.LS[, m] == group.index[g]]
			x.g = x[z.samples.LS[, m] == group.index[g], ]
			r.g = r.samples.LS[(group.index[g]+1),,m]
			q.g = sum(r.g)
			omega.g = omega.samples[z.samples.LS[, m] == group.index[g], m]
			if(q.g >0){
				x.g.1 = x.g[, r.g==1]
				Omega.var.g.inv = diag(1/omega.g)
				Lambda.g = t(x.g.1)%*%Omega.var.g.inv%*%x.g.1				
				if(any(class(tryCatch(solve(Lambda.g), error = function(e) e)) == "error"))
					Lambda.g.inv = solve(Lambda.g + 0.01 * diag(ncol(Lambda.g)))
				else
					Lambda.g.inv = solve(Lambda.g)

	            beta.est.tmp = Lambda.g.inv%*%t(x.g.1)%*%Omega.var.g.inv%*%cbind(y.g)
	           	beta.est[r.g==1, ] = beta.est.tmp
	           	beta.est.group[m, , g] = beta.est
	           	est.mean.g = x.g.1%*%beta.est.tmp
	            SSE = t(y.g - est.mean.g)%*%Omega.var.g.inv%*% (y.g - est.mean.g) 
	        }
	        else{
	        	SSE = (sum((y.g)^2))	
	        }
	        
	        sigma.est.group[m, g] = est.std = sqrt(SSE/ length(y.g)) 
			for(i in 1:num.of.obs){
				est.std.i = sqrt(omega.samples[i, m])*est.std
				est.mean.i = x[i, , drop=FALSE] %*% beta.est
				prob.cube[m, i, g] = rho.samples.LS[(group.index[g]+1),m]*dnorm(y[i], est.mean.i, est.std.i)
				prob.cube.ICL[m, i, g] = dnorm(y[i], est.mean.i, est.std.i)
			}	
		}
		
		prob.cube.cluster[m,,] = t(apply(prob.cube[m,,], 1, function(x) x/sum(x))) 
	}

loglike = 0
loglike.ICL = 0
for(m in 1:num.of.iter.eff)
	for(sub.index in 1:num.of.obs){
		loglike = loglike + log(sum(prob.cube[m, sub.index, ]))
		loglike.ICL = loglike.ICL + log((prob.cube.ICL[m, sub.index, which.max (prob.cube.cluster[m, sub.index, ])]))
	}
AIC = -2*loglike/num.of.iter.eff + 2*num.of.obs
BIC = -2*loglike/num.of.iter.eff + log(num.of.obs)*num.of.groups*(2 + num.of.covariates)
ICL = -2*loglike.ICL/num.of.iter.eff + log(num.of.obs)*num.of.groups*(2 + num.of.covariates)


library(label.switching)
r.samples.rm.nan = array(0, c(num.of.groups, num.of.covariates, num.of.iter.eff))
z.samples.LS.after = z.samples.LS
permutation.matrix = stephens(prob.cube)$permutations

for(m in 1:num.of.iter.eff){
	r.samples.rm.nan[, , m] = r.samples.LS[,,m][complete.cases(r.samples.LS[,,m]*0), ,drop=FALSE]
	r.samples.rm.nan[, , m] = r.samples.rm.nan[permutation.matrix[m,], , m]

	group.index = unique(z.samples.LS[, m])
	#if(is.unsorted(permutation.matrix[m,]))
		for(g in 1:num.of.groups)
			z.samples.LS.after[z.samples.LS[, m]==group.index[g], m] = g 
	#else
		#z.samples.LS.after[, m] = 1:num.of.groups
}
group.level = 0:(num.of.groups-1)
z.samples.LS.after.Relabel = z.samples.LS.after
relabel.times = 0 
for(m in 1:num.of.iter.eff)
{
	if(is.unsorted(permutation.matrix[m,], strictly=TRUE)){
		a = which(z.samples.LS.after[, m] == 0)
		b = which(z.samples.LS.after[, m] == 1)
		z.samples.LS.after.Relabel[a] = 1
		z.samples.LS.after.Relabel[b] = 0
		relabel.times = relabel.times + 1
	}

}




z.est = unlist(lapply(apply(z.samples.LS.after.Relabel, 1, table),  function(x) names(which.max(x))))


sigma.est.group.relabel = sigma.est.group
for(i in 1:nrow(permutation.matrix))
	sigma.est.group.relabel[i, ] = sigma.est.group[i, permutation.matrix[i,]]


beta.est.group.LS.after = beta.est.group
for(j in 1:num.of.covariates)
	for(perm.index in 1:nrow(permutation.matrix))
		beta.est.group.LS.after[perm.index, j, ] = beta.est.group[perm.index, j, permutation.matrix[perm.index,]]



if(0){
beta.est.credit = apply(beta.est.group.LS.after, c(2, 3), mean)
rownames(beta.est.credit) = colnames(x)
#pdf(file="CreditGrowth_beta_normal_1021.pdf", height = 14, width = 7)
par(las=2)
par(mar=c(5,14,4,2))

barplot(t(beta.est.credit), main= expression(paste("Posterior estimate of ", beta)), xlim=c(min(beta.est.credit)-2, max(beta.est.credit)+2), 
		xlab=expression(beta), col=2:(num.of.groups+1),
		names.arg = rownames(beta.est.credit), beside=TRUE, horiz=TRUE)

dev.off()


par(mfrow=c(2, 1))
plot(beta.est.group[, 9, 1], type="l", ylim = range(beta.est.group[, 9, ]), col=2)
lines(beta.est.group[, 9, 2], col=3)

plot(beta.est.group.LS.after[, 9, 1], type="l", ylim = range(beta.est.group[, 9, ]), col=2)
lines(beta.est.group.LS.after[, 9, 2], col=3)
}


#=======================================================================#

#=======================================================================#

m.keep = numeric(0)
for(m in (num.of.iter.eff.remove+1):num.of.iter.eff){
	r.tmp = r.samples.LS[,,m][complete.cases(r.samples.LS[,,m]*0), ,drop=FALSE]
	if(nrow(r.tmp) == num.of.groups ){
		m.keep = c(m.keep, m)
		r.samples.rm.nan[, , m-num.of.iter.eff.remove] = r.samples.LS[,,m][complete.cases(r.samples.LS[,,m]*0), ,drop=FALSE]
		rho.samples.rm.nan[, m-num.of.iter.eff.remove] = rho.samples.LS[,m][complete.cases(r.samples.LS[,,m]*0), drop=FALSE]
		#similarity = apply(r.samples.rm.nan[, , m-num.of.iter.eff.remove], 1, function(x) sum(x^2))
		#if( sum( (r.samples.rm.nan[2, , m-num.of.iter.eff.remove] - r.samples.rm.nan[1, , num.of.iter.eff.remove])^2 )  < sum( (r.samples.rm.nan[1, , m-num.of.iter.eff.remove] - r.samples.rm.nan[1, , num.of.iter.eff.remove])^2 ) )
			r.samples.rm.nan[, , m-num.of.iter.eff.remove] = r.samples.rm.nan[order(rho.samples.rm.nan[, m-num.of.iter.eff.remove]), , m-num.of.iter.eff.remove] 
		#if(num.of.obs.in.groups[1]<num.of.obs.in.groups[2])
			#r.samples.rm.nan[, , m-num.of.iter.eff.remove] = r.samples.rm.nan[c(2,1), , m-num.of.iter.eff.remove]
		
	}
}

r.samples.intercept = r.samples.rm.nan

#for(i in 1:num.of.iter.eff)
#	r.samples.intercept[,,i] = r.samples.intercept[order(beta.est.group.LS.after)[i, 1, ], ,i]


r.est = apply(r.samples.intercept, c(1, 2), mean)
colnames(r.est) = colnames(x)



PIP.Beta.all = cbind(t(r.est), apply(beta.est.group.LS.after, c(2, 3), mean))

r.est.dummy = r.est 
r.est.dummy[r.est<0.5] = 0
r.est.dummy[r.est>0.5] = 1


r.est = as.table(r.est)


r.est.simple = r.est

r.est.simple[r.est<0.5] = 0

colnames(r.est.simple)[apply(r.est.simple, 2, function(x) all(x==0))] = ""

r.est.simple.show = r.est.simple[,!apply(r.est.simple, 2, function(x) all(x==0))]


beta.est.credit = apply(beta.est.group.LS.after, c(2, 3), mean)
rownames(beta.est.credit) = colnames(x)

beta.est.credit[t(r.est)<0.5] = 0

beta.est.credit[!apply(r.est.simple, 2, function(x) all(x==0)), ]

r.beta.est.simple = cbind(t(r.est.simple[,!apply(r.est.simple, 2, function(x) all(x==0))]), beta.est.credit[!apply(r.est.simple, 2, function(x) all(x==0)), ])

print(xtable((r.beta.est.simple),floating=FALSE,latex.environments=NULL,booktabs=TRUE), file = "/Users/kjlee/Research/KJLEE_Papers/MixtureRegressionModels/YC_Martin_KJ/r.tex")




par(las=2)
par(mar=c(5,14,4,2))

barplot(r.est, main= expression(paste("Posterior estimate of P(", gamma, "=1)")), 
xlab="Probability", col=2:(num.of.groups+1),
	names.arg = colnames(r.est), beside=TRUE, horiz=TRUE)

abline(v=0.5, lwd=2, col=2)


r.est = apply(r.samples.rm.nan, c(1, 2), mean)
colnames(r.est) = colnames(x)

r.est.dummy = r.est 
r.est.dummy[r.est<0.5] = 0
r.est.dummy[r.est>0.5] = 1


r.est = as.table(r.est)

par(las=2)
par(mar=c(5,14,4,2))

barplot(r.est, main= expression(paste("Posterior estimate of P(", gamma, "=1)")), 
xlab="Probability", col=2:(num.of.groups+1),
	names.arg = colnames(r.est), beside=TRUE, horiz=TRUE)


pdf(file="/Users/kjlee/Research/KJLEE_Papers/MixtureRegressionModels/YC_Martin_KJ/Figures/CreditGrowth_gamma_robust_RM.pdf", height = 18, width = 7)

par(las=2)
par(mar=c(5,14,4,2))


r.est.simple = r.est

r.est.simple[r.est<0.5] = 0

colnames(r.est.simple)[apply(r.est.simple, 2, function(x) all(x==0))] = ""

barplot(r.est.simple, main= expression(paste("Posterior estimate of P(", gamma, "=1)")), 
xlab="Probability", col=2:(num.of.groups+1),
	names.arg = colnames(r.est.simple), beside=TRUE, horiz=TRUE)


#barplot(r.est, main= expression(paste("Posterior estimate of P(", gamma, "=1)")), 
#  xlab="Probability", col=c("red"), names.arg = colnames(r.est), beside=TRUE, horiz=TRUE)

abline(v=0.5, lwd=2, col=2)

dev.off()


beta.est.credit = apply(beta.est.group.LS.after, c(2, 3), mean)
rownames(beta.est.credit) = colnames(x)

beta.est.credit[t(r.est)<0.5] = 0
pdf(file="/Users/kjlee/Research/KJLEE_Papers/MixtureRegressionModels/YC_Martin_KJ/Figures/CreditGrowth_beta_robust_RM.pdf", height = 18, width = 7)
par(las=2)
par(mar=c(5,14,4,2))

barplot(t(beta.est.credit), main= expression(paste("Posterior estimate of ", beta)), xlim=c(min(beta.est.credit)-2, max(beta.est.credit)+2), 
		xlab=expression(beta), col=2:(num.of.groups+1),
		names.arg = colnames(r.est.simple), beside=TRUE, horiz=TRUE)

dev.off()

#pdf(file="CreditGrowth_gamma_normal_1021.pdf", height = 14, width = 7)


pdf(file="Credit_Omega_barplot_1110.pdf", height = 14, width = 7)
par(las=2)
par(mar=c(5,14,4,2))

barplot(apply(omega.samples, 1, mean) , main= expression(paste("Posterior estimate of ", omega)), 
xlab=expression(omega), names.arg = rownames(X), horiz=TRUE)
abline(v=3, lwd=2, col=2)

dev.off()

pdf(file="Credit_Omega_barplot_thin.pdf", height = 14, width = 7)
par(las=2)
par(mar=c(5,14,4,2))

barplot(apply((omega.samples[, seq(1, length(omega.samples[1, ]), by=100)]), 1, mean), main= expression(paste("Posterior estimate of ", omega)), 
xlab=expression(omega), names.arg = rownames(X), horiz=TRUE)
abline(v=3, lwd=2, col=2)

dev.off()


}


if(0)
{
country.of.interest = rownames(X)

#country.of.interest = unlist(lapply(country.of.interest.abb, isoToName))

#rwmGetISO3

gpclibPermit()
world.map <- readOGR(dsn="/Users/kjlee/Research/KJLEE_Papers/MixtureRegressionModels/Data/TM_WORLD_BORDERS_SIMPL-0.3", layer="TM_WORLD_BORDERS_SIMPL-0.3")
world.ggmap <- fortify(world.map, region = "NAME")
#all.country = unique(world.ggmap$id)
all.country = as.vector(world.map$ISO3)
#all.country.abb = unlist(lapply(all.country, rwmGetISO3))

#all.country %in% country.of.interest
#sum(all.country %in% country.of.interest)
#sort(country.of.interest)
#growth.data.match = match(all.country, country.of.interest)
growth.data.match = match(all.country[all.country %in% country.of.interest], country.of.interest)


id.country = unlist(lapply(all.country, isoToName))
id.country[which(id.country == "United States of America")] = "United States"
n <- length(world.map$ISO3)
df <- data.frame(id = id.country, Grouping = rep(1,n))

## noise
df[!(all.country %in% country.of.interest), 2] =  NA
df[(all.country %in% country.of.interest), 2] = z.est[growth.data.match] #apply(posterior.samples$omega.samples, 1, mean)[growth.data.match] #posterior.samples$z[growth.data.match, num.of.iteration]
#apply(posterior.samples$omega.samples, 1, mean)[growth.data.match]
#df[c(sample(1:100,40)),c("growth", "category")] <- NA

#pdf(file="CreditCrowth_CUMLOSS_WorldMap_1107.pdf", width=14)

ggplot(df, aes(map_id = id)) +
     geom_map(aes(fill = Grouping, color = category), map =world.ggmap, colour = 'black') +
     expand_limits(x = world.ggmap$long, y = world.ggmap$lat) + 
     scale_fill_manual(values = c("red", "green"), labels = c("Group 1", "Group 2"))

ggsave("/Users/kjlee/Research/KJLEE_Papers/MixtureRegressionModels/Kuo-Jung/Latex/Figures/CreditCrowth_CUMLOSS_WorldMap_Group_1108.pdf", width=14, height=7)

sPDF = joinCountryData2Map(df, joinCode='NAME' , nameJoinColumn='id' , verbose='TRUE')


dev.off()

}


if(ClusteringMap){
country.of.interest = rownames(X)

#country.of.interest = unlist(lapply(country.of.interest.abb, isoToName))

#rwmGetISO3

gpclibPermit()
world.map <- readOGR(dsn="/Users/kjlee/Research/KJLEE_Papers/MixtureRegressionModels/Data/TM_WORLD_BORDERS_SIMPL-0.3", layer="TM_WORLD_BORDERS_SIMPL-0.3")
world.ggmap <- fortify(world.map, region = "NAME")
#all.country = unique(world.ggmap$id)
all.country = as.vector(world.map$ISO3)
#all.country.abb = unlist(lapply(all.country, rwmGetISO3))

#all.country %in% country.of.interest
#sum(all.country %in% country.of.interest)
#sort(country.of.interest)
#growth.data.match = match(all.country, country.of.interest)
growth.data.match = match(all.country[all.country %in% country.of.interest], country.of.interest)

id.country = unlist(lapply(all.country, isoToName))
id.country[which(id.country == "United States of America")] = "United States"
n <- length(world.map$ISO3)
df <- data.frame(id = id.country, Group = rep(1,n))

## noise
df[!(all.country %in% country.of.interest), 2] =  NA
df[(all.country %in% country.of.interest), 2] = z.est[growth.data.match] #apply(posterior.samples$omega.samples, 1, mean)[growth.data.match] #posterior.samples$z[growth.data.match, num.of.iteration]
#apply(posterior.samples$omega.samples, 1, mean)[growth.data.match]
#df[c(sample(1:100,40)),c("growth", "category")] <- NA

pdf(file="CreditCrowth_CUMLOSS_WorldMap_Group.pdf", width=14)

ggplot(df, aes(map_id = id)) +
     geom_map(aes(fill = Group, color = category), map =world.ggmap, colour = 'black') +
     expand_limits(x = world.ggmap$long, y = world.ggmap$lat) 
#scale_fill_gradient(low = "red", high = "blue", guide = "colorbar", name="Membership")


sPDF = joinCountryData2Map(df, joinCode='NAME' , nameJoinColumn='id' , verbose='TRUE')

dev.off()
}
#print(sum(posterior.samples$z[, 1] == posterior.samples$z[, num.of.iteration])/num.of.obs)
#print(apply(posterior.samples$r.samples, c(1, 2), mean))
#print(apply(posterior.samples$rho.samples, 1, mean))

