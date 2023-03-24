set.seed(12345)
model_file_name<-'base_mod.txt' #model filename, change to temp_abs_mod.txt for the migration one
number<-10 # percentage of the sample that we use here

print(paste("the number is:",number))

#call packages
library(rjags)
library(beepr)
library(sp)
library(raster)

# read data
load('RUGR - BU LAB PROJECT - detections & ABMI covariates.Rdata')
# sort it
sorted <- match(apply(covariate[, c(1, 14)], 1, paste, collapse = ""), apply(detection[, c(1, 2)], 1, paste, collapse = ""))
detesort<-detection[sorted,]
## subsampling and covariates selection
fine<-grep(150,colnames(covariate))
coarse<-grep(564,colnames(covariate))
finec<-colnames(covariate)[fine]
coarsec<-colnames(covariate)[coarse]
#remove comment for audio wanings during script #beep(2)
n_cov<-number
sa_size<-ceiling(length(detection[,1])/100*number) 
print(paste("Sampling",sa_size,"observations out of ",length(detection[,2])))
print(date())
covs<-sample(1:length(finec),n_cov)
smallsa<-sample(1:length(detesort[,3]),sa_size)

detection<-detesort[smallsa,]
detections<-detection[,3]

sites<-unique(detesort[smallsa,1])
site<-0
for(q in 1:sa_size) site<-c(site,which(sites==detection[q,1]))
site<-site[-1]

obspersite<-rep(0,length(sites))
for(q in 1:length(sites)) obspersite[q]<-length(which(detection[,1]==sites[q]))

print(paste("There are from ",min(obspersite),"to",max(obspersite),"sampling events for each of the",length(sites),"sites"))
fine_covs<-covariate[smallsa,fine[covs]]
coarse_covs<-covariate[smallsa,coarse[covs]]

cova<-fine_covs

#coords

coordata<-covariate[smallsa,c(16,15)]

xs<-ys<-rep(0,length(sites))

for(q in 1:length(sites)) 
	{
	xs[q]<-coordata[which(site==q)[1],2]
	ys[q]<-coordata[which(site==q)[1],1]
}

###### covariates z scores

for(q in 1:n_cov) 
	{
	print(paste("covariate #",q,": ",colnames(cova)[q],", has ",length(unique(cova[,q])),sep=""," unique values"))
	appo<-(cova[,q]-mean(cova[,q],na.rm=TRUE))/sd(cova[,q],na.rm=TRUE)
	cova[,q]<-appo
	}

covar=array(dim=c(length(sites),max(obspersite),n_cov))

###### declaring mcmc list

list_a<-list(n_obs=length(detections),n_cov=length(cova[3,]),det = matrix(nrow=length(sites), ncol = max(obspersite),data=0),cov=covar,sites=length(sites),obspersite=obspersite)


for(q in 1:list_a$sites) 
	{
	testu<-detection[which(site==q),]
	teind<-which(site==q)
	for(w in 1:obspersite[q]) list_a$det[q,w]<-testu[w,3]
	for(w in 1:n_cov) covar[q,1:obspersite[q],w]<-cova[teind,w]
	}

covar[is.na(covar)==TRUE]<-0
list_a$cov<-covar
	

phinit<-rep(0,2)

# get some randomly generated initial values

z_init<-gammini<- matrix(NA, nrow = list_a$sites, ncol = max(obspersite))
for (q in 1:list_a$sites) {
  z_init[q, 1:obspersite[q]] <- 1+rpois(obspersite[q],runif(1,0,1))
  gammini[q, 1:obspersite[q]]<-runif(obspersite[q],0,1)
  phinit[q]~runif(1,0,1)
}
inits_a <- list(upperbound=1,alpha = rnorm(1), lambda = runif(1, 0, 1), beta = rnorm(list_a$n_cov), betasite = rnorm(list_a$sites), z = z_init)

#### function here because I need it somewhere
invlogit <- function(x) {
  exp(x) / (1 + exp(x))
}


######## MCMC parameters
iterations<-24e3
thinning<-4
n_chains<-3
tot_iter_per_chain<-iterations/thinning
tot_iter<-tot_iter_per_chain*n_chains
burn<-iterations/2
track <- c("beta","alpha","betasite","lambda","z","upperbound","chi","psi","phi","gamma","tempabs","sitabs")
########### MCMC declaration 
print(paste("declaring the model with",n_chains,"chains"))
print(date())
model_a <- jags.model(file=model_file_name, data = list_a, n.chains = n_chains,inits=inits_a)
#remove comment for audio wanings during script #beep(2)
print("the model has been declared")
print(date())
########### MCMC burn-in 
update(model_a,burn)
print(paste("the first",burn,"iterations have been discarded"))
print(date())
#remove comment for audio wanings during script #beep(6)
print(paste("starting the sampling with",iterations,"iterations retaining one every",thinning))
print(date())
########### MCMC posterior sampling
post_draws<-coda.samples(model=model_a,n.iter=iterations,variable.names=track,thin=thinning)
#remove comment for audio wanings during script #beep(8)
post_summ<-summary(post_draws) 
####### make a structure with those variables entirely positive or negative (i.e. not overlapping zero)
interesting<-post_summ[[2]][which((post_summ[[2]][,1]>0)|(post_summ[[2]][,5]<0)),]
######## remove the zetas as they will always be positive | need to do this with other parameters if not interested
ind<-grep(paste("z\\[",q,sep=""),rownames(interesting))
interesting<-interesting[-ind,]

print(paste("summaries of the posterior distributions for",track,"are in the post_summ tables ([[1]] for summary statistics and [[2]] for quantiles), while the variables entirely above or below zero are in the interesting table"))
print(date())
### get mean abundance per site 

meanAb<-rep(0,list_a$sites)
for(q in 1:list_a$sites)
	{
	ind<-grep(paste("z\\[",q,sep=""),colnames(post_draws[[1]]))
	meanAb[q]<-0
	for(w in 1:n_chains) meanAb[q]<-meanAb[q]+sum(post_draws[[w]][,ind])/(tot_iter_per_chain*length(ind))
	meanAb[q]<-meanAb[q]/n_chains
	}

coords <- data.frame(x =xs, y =ys, meanab=meanAb) 
remove<-unique(c(which(is.na(coords[,1])),which(is.na(coords[,2]))))
coords<-coords[-remove,]


# Convert to spatial points data frame
coordinates(coords) <- ~x+y
proj4string(coords) <- CRS("+init=epsg:26912") # replace with your coordinate system

# Create raster extent
corners <- matrix(nrow=2,ncol=2)
corners[1,1]<-min(coords@coords[,1],na.rm=TRUE)
corners[1,2]<-max(coords@coords[,1],na.rm=TRUE)
corners[2,1]<-min(coords@coords[,2],na.rm=TRUE)
corners[2,2]<-max(coords@coords[,2],na.rm=TRUE)
rastextent <- extent(corners) + 1

# Convert to raster
rasres <- 1
raster_abundance <- raster(rastextent, resolution = rasres)
raaa<-rasterize(coords, raster_abundance, field = "meanab")
	

# Plot the raster
print("plotting mean estimated abundance map")
print(date())
plot(raaa, col = rev(heat.colors(10)), main = "mean estimated abundance map",xlim=c(corners[1,1],corners[1,2]),ylim=c(corners[2,1],corners[2,2]),xlab="long",ylab="lat")

	

#remove comment for audio wanings during script #beep(4)

########### covariate

#which one is alpha
alphind<-grep("alpha",colnames(post_draws[[1]]))
alpha_pm <- 0

for(w in 1:n_chains) alpha_pm<-alpha_pm+sum(post_draws[[w]][,alphind])/(tot_iter_per_chain*length(alphind))
alpha_pm<-alpha_pm/n_chains

#betas
beta_pm<-rep(0,n_cov)
for(q in 1:n_cov)
	{
		ind<-grep(paste("beta\\[",q,"\\]",sep=""),colnames(post_draws[[1]]))
		beta_pm[q]<-0
		for(w in 1:n_chains) beta_pm[q]<-beta_pm[q]+sum(post_draws[[w]][,ind])/(tot_iter_per_chain*length(ind))
		beta_pm[q]<-beta_pm[q]/n_chains
	}

betasite_pm<-rep(0,list_a$sites)
for(q in 1:list_a$sites)
	{
		ind<-grep(paste("betasite\\[",q,"\\]",sep=""),colnames(post_draws[[1]]))
		betasite_pm[q]<-0
		for(w in 1:n_chains) betasite_pm[q]<-betasite_pm[q]+sum(post_draws[[w]][,ind])/(tot_iter_per_chain*length(ind))
		betasite_pm[q]<-betasite_pm[q]/n_chains
	}

################# covariate_map | right now this only works with 10 covariates

cov_number<-5

covas<-1:number

covas<-covas[-cov_number]

post_mean_beta <- beta_pm[cov_number]


p_beta <- matrix(0, nrow = list_a$sites, ncol = max(obspersite))

for (q in 1:list_a$sites) {
  for (w in 1:obspersite[q]) {
    p_beta[q, w] <- invlogit(alpha_pm + betasite_pm[q] + post_mean_beta * list_a$cov[q, w, cov_number] + beta_pm[covas[1]] * list_a$cov[q, w, covas[1]] + beta_pm[covas[2]] * list_a$cov[q, w, covas[2]] + beta_pm[covas[3]] * list_a$cov[q, w, covas[3]] + beta_pm[4] * list_a$cov[q, w, 4] + beta_pm[covas[5]] * list_a$cov[q, w, covas[5]] + beta_pm[covas[6]] * list_a$cov[q, w, covas[6]] + beta_pm[covas[7]] * list_a$cov[q, w, covas[7]] + beta_pm[covas[8]] * list_a$cov[q, w, covas[8]] + beta_pm[covas[9]] * list_a$cov[q, w, covas[9]])
  }
}
flatbeta<-rep(0,list_a$sites)
for(q in 1:list_a$sites) flatbeta[q]<-mean(p_beta[q,])

rasterbe <- raster(rastextent, resolution = rasres)

becoords <- data.frame(x = xs, y = ys, value = flatbeta)
scarti<-0
for(q in 1:3) scarti<-c(scarti,which(is.na(becoords[,q])))
scarti<-scarti[-1]
becoords<-becoords[-scarti,]

# Convert to spatial points data frame
coordinates(becoords) <- ~x+y

crs(rafbe) <- CRS("+init=epsg:26912")
rafbe <- rasterize(becoords, rasterbe, field = "value")

save.image(version=2,file=paste(model_file_name,"_nm_",substr(date(),9,10),substr(date(),5,7),".RData",sep=""))