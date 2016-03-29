#################################################
# This R code corresponds to the manuscript
# to replicate 
#
# ******     Sum of left and right ICV  ***********
#
# results and figures from SIMULATED DATA
# hence these results will not reflect exact replicates from results in manuscript

# Comparisons of neurodegeneration over time between healthy ageing and Alzheimer's 
# disease cohorts
#
# *** insert ePubs or journal reference ***


#################################
# required libraries
library(rjags)
library(coda)
require(ggplot2)
require(lme4)
library(MASS)
library(reshape2)
require(quantregGrowth)

# Load simulated data - generated by Simulate_MRI_data.r
load("final.sim.dat.Rdata")
sim.d<- final.sim.dat

noObs<- dim(sim.d)[1]
noObs
sim.d$ID<- factor(as.character(sim.d$ID))
ppl<- length(levels(sim.d$ID)) # get 240 people
ppl
sim.d$Age_MCI<- sim.d$mci.bin*sim.d$StndAge
sim.d$Age_AD<- sim.d$ad.bin*sim.d$StndAge

meanAge<- 75; sdAge<- 7  # These are to convert back to Age scale - as described in
# Section 2 of manuscript
head(sim.d)

person<- as.numeric(as.character(sim.d$ID))
##############################################################################################
##############################################################################################

############################
#  models with interaction #
############################

#-----------------------------
# vent/SumH~ stndAge + x1(mci) + x2(ad) + stndAge*x1(mci) + stndAge*x2(AD) 
#--------------------------------

data.inter<- list(N = ppl, Y = sim.d$Hippo.ICV, 
                  age=sim.d$StndAge, 
                  x1= sim.d$mci.bin, x2=sim.d$ad.bin, predY=rep(NA, noObs), 
                  person=person, obs=noObs,
                  Age_MCI=sim.d$Age_MCI, 
                  Age_AD= sim.d$Age_AD,
                  #  multivariate random effects model
                  prec = structure(.Data = c(0.0001,0     ,0      ,0,
                                             0     ,0.0001,0      ,0,
                                             0     ,0     , 0.0001,0,
                                             0     ,0     ,0      ,0.0001), .Dim = c(4,4)),
                  mean = c(0,0,0,0), R = structure(.Data = c(0.0001,0 ,0 ,0,
                                                             0 ,0.0001,0 ,0,
                                                             0 ,0 ,0.0001,0,
                                                             0 ,0 ,0 ,0.0001), .Dim = c(4,4)))
str(data.inter)

response<- 'SumHippo'
y.response<- sim.d$Hippo.ICV

###############
## Run model ##
###############

model <- jags.model('BLMEmodel_file.bug', data = data.inter, n.chains = 2, n.adapt = 100)
#model
n.iter <- 200000; thin <- 50; burn.in<- 100000

tme<- proc.time()
update(model, n.iter) 
samp.2<- coda.samples(model, c('sigma2', 'mu.beta', 'omega', 'beta2', 'beta3'), n.iter = n.iter, start = burn.in, thin = thin)
(proc.time() - tme)/60
# took 17 minutes on home computer

summary(samp.2, start = burn.in)
effectiveSize(samp.2) # such low ESS: notw diagnostics look ok

##########################################
## retreive MCMC chains for parameters  ##
##########################################

mcmc <- samp.2
# convert mcmc chains to data.frame object 
for (j in 1:2) {
  mcmc[[j]] <- as.data.frame(mcmc[[j]])
  n <- dim(mcmc[[j]])[1]
  mcmc[[j]][,"id"] <- 1:n
  mcmc[[j]][,"chain"] <- rep(j,n)
}

mcmcs <- rbind(mcmc[[1]], mcmc[[2]])

# Additional: convert covariance matrix onto correlation matrix to interpret the relationship between the random effects
#PostMeans<- colMeans(mcmcs)
#CovarianceM<- matrix(c(PostMeans[7], PostMeans[11], PostMeans[15], PostMeans[19],
#                       PostMeans[8], PostMeans[12], PostMeans[16], PostMeans[20],
#                       PostMeans[9], PostMeans[13], PostMeans[17], PostMeans[21],
#                       PostMeans[10], PostMeans[14], PostMeans[18], PostMeans[22]) , nrow=4, ncol=4)
#CorrM<- cov2cor(CovarianceM)
#CorrM

# note: we only need the bottom diagonal of Omega matrix
# no need to change whole bunch of syntax for all inferences! Just change the format of the mcmcs object
mcmcs<- mcmcs[, -c(11,15,19,16,20,21)]
colnames(mcmcs)<- c("beta2", "beta3",  "beta.0", "beta.1", "beta.4", "beta.5", "omega11", "omega21", "omega31", "omega41", "omega22", "omega32", 
                    "omega42", "omega33", "omega43", "omega44", "sigma2", "id", "chain")

mcmcs.Hippo<- mcmcs
save(mcmcs.Hippo, file="mcmcs.Hippo.Rdata")
#load("mcmcs.Hippo.Rdata")

##_____________________________________________________________________________
##_____________________________________________________________________________
                      ############################################
                      # Analysis 1                               #
                      # multiple comparisons among population    #
                      # groups.                                  #
                      # Section 3.3 and 4.1 of manuscript        #
                      ############################################ 

hc<- density(mcmcs$beta.1)
mci<- density(mcmcs$beta.1 + mcmcs$beta.4)
ad<- density(mcmcs$beta.1 + mcmcs$beta.5)
data.1<- data.frame(x= c(hc$x, mci$x, ad$x), y=c(hc$y, mci$y, ad$y), cls=rep(c("HC","MCI", "AD"), each=512))

#######################################################
# posterior mean and CI for each diagnosis
# Similar results for Table 1
#######################################################

chains.diag<- data.frame(hc.chain = mcmcs$beta.1, mci.chain=  mcmcs$beta.1 + mcmcs$beta.4, ad.chain=mcmcs$beta.1 + mcmcs$beta.5 )
mEans<- data.frame(diag = c("hc", "mci", "ad"), post.means = rep(0, 3), ci.high = rep(0, 3), ci.low = rep(0, 3))
for(i in 1:3){
  mEans$post.means[i] = mean(chains.diag[, i])
  mEans$ci.high[i] = quantile(chains.diag[, i], probs = 0.025)
  mEans$ci.low[i] = quantile(chains.diag[, i], probs = 0.975)
}
mEans

###################################################################################################################
# Aside: difference in population means: note: we only want the simulations for each of the densities plotted in Figure 1
# #################################################################################################################
#hc.slope<- mcmcs$beta.1
#mci.slope<- mcmcs$beta.1 + mcmcs$beta.4
#ad.slope<- mcmcs$beta.1 + mcmcs$beta.5

# visualise the differences between the slopes: hc.slope and mci.slope, hc.slope and ad.slope
#source("multiplot.r")
#windows()
#p.1<- ggplot(data = mcmcs, aes(x = beta.4)) + geom_density(alpha = 0.2, fill = "blue") + ggtitle('SIMULATED Estimated difference between HC and CI groups (beta_4)') +theme_bw()
#p.2<- ggplot(data = mcmcs, aes(x = beta.5)) + geom_density(alpha = 0.2, fill = "red") + ggtitle('SIMULATED Estimated difference between HC and AD groups (beta_5)')+ theme_bw()
#multiplot(p.1, p.2)

#######################################################################################
# credible intervals for parameters - recall 
# on code [beta.0, beta.1, beta.4, beta.5] = [beta0, beta3, beta4, beta5] on paper
# bottom of Table 1 in manuscript
#######################################################################################

no.param<- 17
resuLts<- data.frame(param.name = c("beta2", "beta3", "beta.0", "beta.1", "beta.4", "beta.5", "omega11", "omega21",
                                    "omega31", "omega41", "omega22", "omega32", "omega42", "omega33", "omega43", "omega44", "sigma2"),
                     post.means= rep(0, no.param),up.ci= rep(0, no.param),low.ci= rep(0, no.param) )
for(i in 1:no.param){
  resuLts$post.means[i]<- mean(mcmcs[, i])
  resuLts$up.ci[i]<- quantile(mcmcs[, i], probs = 0.025)
  resuLts$low.ci[i]<- quantile(mcmcs[, i], probs = 0.975)
}
resuLts

#######################################################################
# Figure 2 on paper: (bottom density plots)
# for top of Figure 2 see line ~350
#######################################################################
colnames(data.1)[3]<- "Diagnosis"

windows()   
hippo.fixed<- ggplot(data.1, aes(x=x, y=y, colour=Diagnosis)) + geom_line() + 
  ggtitle(paste("SIMULATED ",response, " posterior density plots of global slopes \n across all diagnosis (for all ages)")) + 
  geom_vline(data = mEans, xintercept=c(mEans[1,2], mEans[2,2], mEans[3,2]), linetype="dotted") + 
  xlab("rate of deterioration (slope)") + ylab("density") + theme_bw() +
  scale_colour_discrete(breaks=c("AD","MCI","HC"))  # this makes the ordering of diagnosis AD -> MCI -> HC and retains their colour

hippo.fixed

#################################################################################

# for the standard error of each of these distributions - look at the summary results above!
mci.less.hc<-ad.less.mci<-0

slope.mci<- mcmcs$beta.1 + mcmcs$beta.4
slope.ad<- mcmcs$beta.1 + mcmcs$beta.5
slope.hc<- mcmcs$beta.1

########################################
# compute expression (6) in manuscript
# whose results are in Table 2
########################################
for(i in 1:dim(mcmcs)[1]){
  if(slope.mci[i]< mEans[1,2]){
    mci.less.hc<- mci.less.hc + 1
  }
  if(slope.ad[i] < mEans[2,2]){
    ad.less.mci<- ad.less.mci + 1
  } 
}

# hippo
prob.1 = mci.less.hc/dim(mcmcs)[1]
prob.1 # p(MCI < AD)

prob.2 = ad.less.mci/dim(mcmcs)[1]
prob.2 # P(AD < MCI)

#_________________________________________________________________________________________
#_________________________________________________________________________________________

                    ############################################
                    # Analysis 2                               #
                    # rank participants by randome effects     #
                    # for HC that is B_1i,                     #
                    # for MCI that is B_1i + B4i               #
                    # for AD that is B_1i + B5i                #
                    #                                          #
                    # Section 3.4 and 4.2 of manuscript        #
                    ############################################ 

# note: All the random effects come from a multivariate normal.
rand.eff.Hippo<- coda.samples(model, c('beta'), n.iter=n.iter, thin=thin)

save(rand.eff.Hippo, file = "rand.eff.Hippo.Rdata")  # took about 10 minutes
#load("rand.eff.Hippo.Rdata")
samp.1<- rand.eff.Hippo

#head(samp.1)
mcmc.1 <- samp.1
for (j in 1:2) {
  mcmc.1[[j]] <- as.data.frame(mcmc.1[[j]])
  n <- dim(mcmc.1[[j]])[1]
  mcmc.1[[j]][,"id"] <- 1:n
  mcmc.1[[j]][,"chain"] <- rep(j,n)
}
# these are the chain results for each prediction value
mcmcs.1 <- rbind(mcmc.1[[1]], mcmc.1[[2]])
#head(mcmcs.1)
dim(mcmcs.1) # we have 260 people

# recall beta = [b0, b1, b4, b5]
beta1.sims<- mcmcs.1[, 261:520]
beta4.sims<- mcmcs.1[, 521:780]
beta5.sims<- mcmcs.1[, 781:1040]

######################################################################################
# recall from Simulated_MRI_data.r, there were 157 HC, 34 MCI and 42 AD
# as well as 27 Converters of various forms (HC -> AD, MCI -> AD)

#load("diag.Rdata")
diag2<- cbind(ID=seq(1:260), diag)

# get baseline labels - but recall this is unbalanced
ind.hc<- c(seq(1:157), seq(from=234, to=245, by=1))
ind.mci<- c(seq(from=158, to=191, by=1), seq(from=246, to=260, by=1))
ind.ad<- seq(from=192, to=233, by=1)

# from the above we see the ID orders should be (to include converters)
id.all<- c(ind.hc, ind.mci, ind.ad)

mci.slop.f<- beta1.sims[,c(ind.mci)] + beta4.sims[, c(ind.mci)]; colnames(mci.slop.f)<- ind.mci
ad.slop.f<- beta1.sims[,c(ind.ad)] + beta5.sims[, c(ind.ad)];colnames(ad.slop.f)<- ind.ad
hc.slop.f<- beta1.sims[, c(ind.hc)]; colnames(hc.slop.f)<- ind.hc

# From Simulate_MRI_data.r recall ordering of participants
part.diag.order<- c(rep("HC", 157), rep("Converter", 12), 
                    rep("MCI", 34), rep("Converter", 15), 
                    rep("AD", 42))

hc.mci.ad.chains<- cbind(hc.slop.f, mci.slop.f, ad.slop.f)

dat.s<- data.frame(fin.slop.mean = colMeans(hc.mci.ad.chains),
                   up.ci= apply(hc.mci.ad.chains, 2, quantile, probs=0.975),
                   low.ci = apply(hc.mci.ad.chains, 2, quantile, probs=0.025),
                   Diag.conv = part.diag.order, ID=id.all)

# rank in order from smallest to largest
dat.s.2<- dat.s[order(dat.s$fin.slop.mean),]

#################################################################################################
# correspond with Table 3 in manuscript - list top and bottom 5 participants by ranked order
#################################################################################################

head(dat.s.2, 5)
tail(dat.s.2, 5)

########################################################
# prep data for boxplots similar to manuscript Figure 4
########################################################

ppl.chains<- hc.mci.ad.chains
chains.ord<- ppl.chains[, as.character(dat.s.2$ID)]

dat.l<- melt(chains.ord)
head(dat.l)

colnames(dat.l)<- c("participant", "rate.of.deterioration")
dat.l$Diagnosis<- rep(dat.s.2$Diag.conv, each = 8000)
to.plot.dat<- dat.l

# recall 8 converters of interest were HC -> MCI whose ID's 234 -> 241
conv.int<- subset(to.plot.dat, participant %in% seq(from=234, to=241, by=1))

how.many.conv<- length(seq(from=234, to=241, by=1))
c.means<- c.min<- c.max<- rep(0, how.many.conv)

for(j in 1:how.many.conv){
  t<- j*dim(ppl.chains)[1]
  c.means[j]<- mean(conv.int[1:t, 2])
  c.min[j]<- min(conv.int[1:t, 2])
  c.max[j]<- max(conv.int[1:t,2])
} 
very.bottomHippo<- rep(c(-0.2, -0.22), 4)  #rep(-0.02, 8)

conv.dat1<- data.frame(participant = factor(seq(from=234, to=241, by=1)), rate.of.deterioration = very.bottomHippo, 
                       Diagnosis = rep("Converter", how.many.conv))

conv.dat2<- data.frame(participant = factor(seq(from=234, to=241, by=1)), rate.of.deterioration = c.means, Diagnosis = rep("Converter", 8), 
                       ymin = very.bottomHippo, ymax = c.means) 


#########################################
# Plot similar to Figure 4 of manuscript
#########################################

windows()
p.inf2<- ggplot(to.plot.dat, aes(x = participant, y = rate.of.deterioration, fill=Diagnosis)) + geom_boxplot(outlier.size = 1) + 
  ggtitle(paste("SIMULATED ",response, " model with interactions \n random slopes in order of increasing deterioration")) + theme_bw()#+
p.inf2 # plot ranked box-plots

# include labels for participants of interest
p.inf2+ geom_text(data=conv.dat1, aes(label = participant), hjust = 1, vjust = 0, fontface = "bold", size = 4, angle = 90, colour = "red") +
  geom_pointrange(data = conv.dat2, aes(x = participant,y = rate.of.deterioration, ymin = ymin, ymax = ymax, colour = "red"),
                  show.legend = FALSE) 


###############################################################################################
# Posterior rank distribution for converters
# note: This analysis can be repeated for a subset of the Simulated data - 
# for example the first three time points - as done in the manuscript
###############################################################################################

head(chains.ord,0)  # this has ordered the MCMC columns by posterior means 
nameS<- colnames(chains.ord,0)# whose respective ID's are stored here 

rank.dist<- matrix(0, nrow = dim(chains.ord)[1], ncol=dim(chains.ord)[2])
#temp.mat<- data.frame(rank = rep(0, 260))
temp.vec<- c()
for(i in 1:dim(chains.ord)[1]){
  temp.vec<- t(chains.ord[i,])
  rank.dist[i,]<- order(temp.vec) # as I am only using row names for ranking
}
class(rank.dist)<- "numeric"

dat.f.rankDist<- data.frame(rank.dist)
colnames(dat.f.rankDist)<- nameS

# recall ID's 234 - 260 were converters in our simulated data
rank.dist.conv<- dat.f.rankDist[, as.character(seq(from=234, to=260, by=1))]

to.plot.dat2<- melt(rank.dist.conv)
windows()
ggplot(data = to.plot.dat2, aes(x = value)) + geom_density() +facet_wrap(facets =~variable, ncol = 6, scales="free") + theme_bw() + 
  xlab("Posterior ranking density plots for all (27) converters") + 
  ggtitle("SIMULATED Hippo.ICV: converter rankings based on estimated rate of volume change (individual estimates)")

#_______________________________________________________________________________________
#_______________________________________________________________________________________

                      ############################################
                      # Analysis 3                               #
                      # Compute probability trajectories for     #
                      # four colume quantiles                    #
                      # Section 3.5 and 4.3 of manuscript        #
                      ############################################ 

tauss2<- c(0.15, 0.25, 0.5, 0.75)

# StndAge             Age
# -1.918367375        60.07939
# -1.26530612         65.01615
# -0.61224490         69.95291
# 0.04081633          74.88966
# 0.77551020          80.44351
# 1.34693878          84.76318

#############################
# Figure 1 in manuscript    #
#############################

windows()  # y=TRUE puts in the data points
p2<- gcrq(Hippo.ICV~ps(StndAge, mon=0, lambda=100, pdiff=2), tau=tauss2, data=sim.d)
plot(p2, pch=20,legend=TRUE, xlim = range(sim.d$StndAge), ylim=range(sim.d$Hippo.ICV), xlab = "", 
     ylab="", lwd=0.5) # this will do the basic plot
par(new=TRUE)
plot(sim.d$StndAge, sim.d$Hippo.ICV, ylab = paste(response, ".ICV"), xlab="StndAge", col=c("blue","seagreen", "red"),   #temp2$cls, 
     xlim = range(sim.d$StndAge), ylim=range(sim.d$Hippo.ICV), cex=0.5, pch=19, main = paste("SIMULATED ",response, " volume vs StndAge"))
abline(v = c(-1.91, -1.27, -0.62, 0.041, 0.776, 1.35), col="gray60", lty=2) # for the 6 age groups, see tbl above
abline(h = max(y.response), col="gray60", lty=2)

# - note: these probabilities must all add to 1 - as you have to be in one of the three states
# these are taken from National rate of dementia 2012 report AIHW2012 (AD rates)
# whose predecessor is # http://www.alzheimersanddementia.com/article/S1552-5260(11)00064-1/pdf + Dementia Australia 2007
# Ward et. al 2012 paper Mild cognitive impairment: disparity of incidence and prevalence estimates (MCI)
# also taken from (ABS) http://www.abs.gov.au/AUSSTATS/abs@.nsf/Lookup/4102.0Main+Features50Dec+2012#10 

# further reference from 1995 - a bit old and won't be included
# http://www.ncbi.nlm.nih.gov/pubmed/7715060
lit.p<- data.frame(Age.Diag=c("p.hc.60", "p.mci.60", "p.ad.60", "p.hc.65", "p.mci.65", "p.ad.65",
                              "p.hc.70", "p.mci.70", "p.ad.70", "p.hc.75", "p.mci.75", "p.ad.75", 
                              "p.hc.80", "p.mci.80", "p.ad.80", "p.hc.85", "p.mci.85", "p.ad.85"), 
                   prob=c(0.945, 0.037,0.018, 
                          0.917, 0.055,0.028,
                          0.859, 0.096,  0.045,
                          0.592,0.333,0.075,
                          0.518, 0.357,0.125,
                          0.496, 0.301,0.203))

hippo.ranges<- data.frame(Age = c("60.top25", "60.50.75", "60.25.50", "60.15.25",
                                  "65.top25", "65.50.75", "65.25.50", "65.15.25",
                                  "70.top25", "70.50.75", "70.25.50", "70.15.25",
                                  "75.top25", "75.50.75", "75.25.50", "75.15.25",
                                  "80.top25", "80.50.75", "80.25.50", "80.15.25",
                                  "85.top25", "85.50.75", "85.25.50", "85.15.25"), 
                          lower=c(0.4318455, 0.418423, 0.04247, 0.39054,
                                  0.4302525, 0.4150818, 0.4000778, 0.3889878,
                                  0.4273551, 0.4102399, 0.3911624, 0.3795268,
                                  0.4230760, 0.4036956, 0.3757708, 0.3623180,
                                  0.4174686, 0.3955275, 0.3541077, 0.3375948,
                                  0.4102640, 0.3853550, 0.3250936, 0.3041123),
                          upper=c(0.4970, 0.438455, 0.418423, 0.40247,
                                  0.4970, 0.4302525, 0.4150818, 0.4000778,
                                  0.4970, 0.4276551, 0.4102399, 0.3911624,
                                  0.4970, 0.4230760, 0.4036956, 0.3757708, 
                                  0.4970, 0.4230760, 0.4036956, 0.3757708,
                                  0.4970, 0.4102640, 0.3853550, 0.3250936) )


#########################################################################################################
# To proceed - need predicted values over ages listed above   - for fixed effects trajectories over age
##########################################################################################################
##############################
# (*) for HC: x1, x2 == 0   
# these sections (*) are ONLY relevant for model vol~ stndAge + cls + stndAge*cls
##############################

mod.hc<- function(age, b0, b1, b4, b5,s0, s21, s31, s41, s1, s32, s42, s4, s43,s5, s){
  mu<- c(b0, b1, b4, b5)
  Sigma<- matrix(data=c(s0, s21, s31, s41,s21, s1, s32, s42,s31, s32, s4, s43,s41, s42, s43, s5), nrow=4, ncol=4)
  beta<- mvrnorm(n=1,mu, Sigma) 
  mew<- beta[1] + beta[2]*age
  y.p<- rnorm(1, mew, s)
  return(y.p)
}

###############################
#  (*) for MCI: x1 = 1, x2 = 0
##############################

mod.mci<- function(age, b0, b1, b4, b5,s0, s21, s31, s41, s1, s32, s42, s4, s43,s5, s, b2){
  mu<- c(b0, b1, b4, b5)
  Sigma<- matrix(data=c(s0, s21, s31, s41,s21, s1, s32, s42,s31, s32, s4, s43,s41, s42, s43, s5), nrow=4, ncol=4)
  beta<- mvrnorm(n=1,mu, Sigma) 
  mew<- beta[1] + b2 + (beta[2] + beta[3])*age
  y.p<- rnorm(1, mew, s)
  return(y.p)
}

###############################
#  (*) for AD: x1 = 1, x2 = 1
###############################

mod.ad<- function(age, b0, b1, b4, b5,s0, s21, s31, s41, s1, s32, s42, s4, s43,s5, s, b3){
  mu<- c(b0, b1, b4, b5)
  Sigma<- matrix(data=c(s0, s21, s31, s41,s21, s1, s32, s42,s31, s32, s4, s43,s41, s42, s43, s5), nrow=4, ncol=4)
  beta<- mvrnorm(n=1,mu, Sigma) 
  mew<- beta[1] + b3 + (beta[2] + beta[4])*age
  y.p<- rnorm(1, mew, s)
  return(y.p)
}

##################################
# posterior predictive checks    #  (this creates the mean predictions across age for all diagnosis)
##################################   (this requires functions above)

agePred<- c(-1.91, -1.27, -0.62, 0.041, 0.776, 1.35)
lengthA<- length(agePred)

load("mcmcs.Hippo.Rdata")
mcmcs<- mcmcs.Hippo

sim.HC<- matrix(data=rep(0, dim(mcmcs)[1]*lengthA), nrow=dim(mcmcs)[1] ,ncol=lengthA)
sim.MCI<- matrix(data=rep(0, dim(mcmcs)[1]*lengthA), nrow=dim(mcmcs)[1] ,ncol=lengthA)
sim.AD<- matrix(data=rep(0, dim(mcmcs)[1]*lengthA), nrow=dim(mcmcs)[1] ,ncol=lengthA)

# go thru each estimate of the joint posterior, and compute agePred estimates
for(k in 1:dim(mcmcs)[1]){
  
  b1<- mcmcs[k, 4]; b0<- mcmcs[k, 3]; s<- mcmcs[k, 17]; s0<- mcmcs[k, 7]; s1<- mcmcs[k,11]
  b2<- mcmcs[k, 1]; b3<- mcmcs[k, 2]; b4<- mcmcs[k, 5];  s4<- mcmcs[k, 14]; b5<- mcmcs[k, 6]; s5<- mcmcs[k, 16];
  # off diagonal covariance parameters
  s21<- mcmcs[k, 8]; s31<- mcmcs[k, 9]; s32<- mcmcs[k,12]; s41<- mcmcs[k,10]; s42<- mcmcs[k, 13]; s43<- mcmcs[k,15];
  
  for(p in 1:lengthA){
    sim.HC[k,p]<- mod.hc(agePred[p],b0, b1, b4, b5,s0, s21, s31, s41, s1, s32, s42, s4, s43,s5, s)
    sim.MCI[k,p]<- mod.mci(agePred[p],b0, b1, b4, b5,s0, s21, s31, s41, s1, s32, s42, s4, s43,s5, s, b2)
    sim.AD[k,p]<- mod.ad(agePred[p], b0, b1, b4, b5,s0, s21, s31, s41, s1, s32, s42, s4, s43,s5, s, b3) 
  }
}

head(sim.HC) # These are posterior estimates 
head(sim.MCI)
head(sim.AD)  

##################################################################
# the following is needed for posterior predictions across ages for: 
# Figure 1 in document
# Figures 4 and 5 (inference 3)

#########
# HC
r<- colMeans(sim.HC)
upper.ci.hc<- rep(0, length(agePred)); lower.ci.hc<- rep(0, length(agePred))

for(q in 1:length(agePred)){
  upper.ci.hc[q]<- quantile(sim.HC[, q], probs=0.025)
  lower.ci.hc[q]<- quantile(sim.HC[, q], probs=0.975)
}
hc.sim<-data.frame(agePred=agePred, obs.m=r, lower.ci=upper.ci.hc,upper.ci=lower.ci.hc) 

##########
# MCI
r.mci<- colMeans(sim.MCI)
upper.ci.mci<- rep(0, length(agePred)); lower.ci.mci<- rep(0, length(agePred))

for(q in 1:length(agePred)){
  upper.ci.mci[q]<- quantile(sim.MCI[, q], probs=0.975)
  lower.ci.mci[q]<- quantile(sim.MCI[, q], probs=0.025)
}
mci.sim<-data.frame(agePred=agePred, obs.m=r.mci, lower.ci=lower.ci.mci, upper.ci=upper.ci.mci) 

######    
# AD
r.ad<- colMeans(sim.AD)
upper.ci.ad<- rep(0, length(agePred)); lower.ci.ad<- rep(0, length(agePred))

for(q in 1:length(agePred)){
  upper.ci.ad[q]<- quantile(sim.AD[, q], probs=0.975)
  lower.ci.ad[q]<- quantile(sim.AD[, q], probs=0.025)
}
ad.sim<-data.frame(agePred=agePred, obs.m=r.ad, lower.ci=lower.ci.ad, upper.ci=upper.ci.ad) 

######################################################################################

bee<- dim(sim.HC)[1]

a.60.chain<- matrix(data=rep(0,bee*3), nrow=bee, ncol=3 )
a.60.chain[, 1]<- sim.HC[, 1]; a.60.chain[, 2]<- sim.MCI[, 1]; a.60.chain[, 3]<- sim.AD[, 1]

a.65.chain<- matrix(data=rep(0,bee*3), nrow=bee, ncol=3 )
a.65.chain[, 1]<- sim.HC[, 2]; a.65.chain[, 2]<- sim.MCI[, 2]; a.65.chain[, 3]<- sim.AD[, 2]

a.70.chain<- matrix(data=rep(0,bee*3), nrow=bee, ncol=3 )
a.70.chain[, 1]<- sim.HC[, 3]; a.70.chain[, 2]<- sim.MCI[, 3]; a.70.chain[, 3]<- sim.AD[, 3]

a.75.chain<- matrix(data=rep(0,bee*3), nrow=bee, ncol=3 )
a.75.chain[, 1]<- sim.HC[, 4]; a.75.chain[, 2]<- sim.MCI[, 4]; a.75.chain[, 3]<- sim.AD[, 4]

a.80.chain<- matrix(data=rep(0,bee*3), nrow=bee, ncol=3 )
a.80.chain[, 1]<- sim.HC[, 5]; a.80.chain[, 2]<- sim.MCI[, 5]; a.80.chain[, 3]<- sim.AD[, 5]

a.85.chain<- matrix(data=rep(0,bee*3), nrow=bee, ncol=3 )
a.85.chain[, 1]<- sim.HC[, 6]; a.85.chain[, 2]<- sim.MCI[, 6]; a.85.chain[, 3]<- sim.AD[, 6]


# combine above
Age.chains<-cbind(a.60.chain,a.65.chain, a.70.chain, a.75.chain, a.80.chain, a.85.chain)
colnames(Age.chains)<- c('hc.60', 'mci.60', 'ad.60','hc.65', 'mci.65', 'ad.65','hc.70', 'mci.70', 'ad.70',
                         'hc.75', 'mci.75', 'ad.75','hc.80', 'mci.80', 'ad.80','hc.85', 'mci.85', 'ad.85')

no.chains<- dim(chains.ord)[1]
Age.c<- data.frame(Age.chains)
Age.chains2<- data.frame(chain.60 = c(Age.c$hc.60, Age.c$mci.60, Age.c$ad.60),
                         chain.65 = c(Age.c$hc.65, Age.c$mci.65, Age.c$ad.65),
                         chain.70 = c(Age.c$hc.70, Age.c$mci.70, Age.c$ad.70),
                         chain.75 = c(Age.c$hc.75, Age.c$mci.75, Age.c$ad.75),
                         chain.80 = c(Age.c$hc.80, Age.c$mci.80, Age.c$ad.80),
                         chain.85 = c(Age.c$hc.85, Age.c$mci.85, Age.c$ad.85),
                         Diag = c(rep("HC", no.chains), rep("MCI", no.chains), rep("AD", no.chains)))
head(Age.chains2)

# ******************** top.25th quantile
top.q<- hippo.ranges[c(1,5,9, 13, 17,21) ,]
rownames(top.q)<- seq(1:6)
top.q
quant<- "75 - 100th Quantile"

#----------
# Compute from the model: P(y_tilde| diagnosis, age)
temp<- data.frame(matrix(data=rep(0, no.chains*6*3), nrow = no.chains*3, ncol=6))

for(j in 1:6){  
  for(u in 1:no.chains*3){   
    if(top.q[j,2] < Age.chains2[u,j] & Age.chains2[u,j] < top.q[j,3]){
      temp[u,j]<- 1} 
  }
}
#head(temp)
#colSums(temp)
temp$Diag<- c(rep("HC", no.chains), rep("MCI", no.chains), rep("AD", no.chains))

# now see what each diagnosis comes as
hc.p<- subset(temp, Diag == "HC")
all.ages.hc.p<- colSums(hc.p[, 1:6])/no.chains

mci.p<- subset(temp, Diag == "MCI")
all.ages.mci.p<- colSums(mci.p[1:6])/no.chains

ad.p<- subset(temp, Diag == "AD")
all.ages.ad.p<- colSums(ad.p[1:6])/no.chains

# put it all into one data frame
p.ytilde.Diag.age<- data.frame(probs.mdl = c(all.ages.hc.p, all.ages.mci.p, all.ages.ad.p), ages = rep( c(60,65,70,75,80,85), 3),
                               diagnosis = c(rep("HC", 6), rep("MCI", 6), rep("AD", 6) ) )   
p.ytilde.Diag.age
p.ytilde<- p.ytilde.Diag.age[order(p.ytilde.Diag.age$ages),]
p.ytilde

#---- now to include probabilities from literature, expression for denominator

denom.again<- c(0,0,0,0,0,0)
d<- c(3,6,9,12,15,18)
for(i in 1:6){
  teMp<- d[i]
  denom.again[i] <-  p.ytilde[teMp-2,1]*lit.p[teMp-2,2] +  p.ytilde[teMp-1,1]*lit.p[teMp-1,2] + p.ytilde[teMp,1]*lit.p[teMp,2] 
}
denom.again 

#----- now estimate final probs

temP<- c( rep(denom.again[1], 3), rep(denom.again[2], 3), rep(denom.again[3], 3), rep(denom.again[4], 3), rep(denom.again[5], 3),
          rep(denom.again[6], 3) )

tEmp<- rep(0, 18)
for(a in 1:18){
  tEmp[a]<- p.ytilde[a, 1]*lit.p[a,2]/temP[a]
}

tEmp<- data.frame(final.p=tEmp, age = p.ytilde$ages, diagnosis=p.ytilde$diagnosis)
tEmp
# ---- plot

tEmp[is.na(tEmp)]<- 2
tEmp

####################################################
# Plots similar to Figure 5 for each quantile 

windows()
ggplot(data=tEmp, aes(x = age, y=final.p, colour=diagnosis)) + geom_point() + geom_line() +
  ggtitle(paste("SIMULATED ", quant)) + ylab("Probability") + xlab("Age") +theme_bw() + theme(legend.position = "none")
# Contact primary author for details on Monte Carlo error intervals
#------------------------------------------------------------------------------------------------------------------------

# ************************    second quantile 50.75 quantile
top.q<- hippo.ranges[c(2,6,10,14, 18,22),]
rownames(top.q)<- seq(1:6)
top.q
quant<- "50 - 75th Quantile"

#----------
# Compute from the model: P(y_tilde| diagnosis, age)
temp<- data.frame(matrix(data=rep(0, no.chains*6*3), nrow = no.chains*3, ncol=6))

for(j in 1:6){  
  for(u in 1:no.chains*3){   
    if(top.q[j,2] < Age.chains2[u,j] & Age.chains2[u,j] < top.q[j,3]){
      temp[u,j]<- 1} 
  }
}
#head(temp)
#colSums(temp)
temp$Diag<- c(rep("HC", no.chains), rep("MCI", no.chains), rep("AD", no.chains))

# now see what each diagnosis comes as
hc.p<- subset(temp, Diag == "HC")
all.ages.hc.p<- colSums(hc.p[, 1:6])/no.chains

mci.p<- subset(temp, Diag == "MCI")
all.ages.mci.p<- colSums(mci.p[1:6])/no.chains

ad.p<- subset(temp, Diag == "AD")
all.ages.ad.p<- colSums(ad.p[1:6])/no.chains

# put it all into one data frame
p.ytilde.Diag.age<- data.frame(probs.mdl = c(all.ages.hc.p, all.ages.mci.p, all.ages.ad.p), ages = rep( c(60,65,70,75,80,85), 3),
                               diagnosis = c(rep("HC", 6), rep("MCI", 6), rep("AD", 6) ) )   
p.ytilde.Diag.age
p.ytilde<- p.ytilde.Diag.age[order(p.ytilde.Diag.age$ages),]
p.ytilde

#---- now to include probabilities from literature, expression for denominator

denom.again<- c(0,0,0,0,0,0)
d<- c(3,6,9,12,15,18)
for(i in 1:6){
  teMp<- d[i]
  denom.again[i] <-  p.ytilde[teMp-2,1]*lit.p[teMp-2,2] +  p.ytilde[teMp-1,1]*lit.p[teMp-1,2] + p.ytilde[teMp,1]*lit.p[teMp,2] 
}
denom.again 

#----- now estimate final probs

temP<- c( rep(denom.again[1], 3), rep(denom.again[2], 3), rep(denom.again[3], 3), rep(denom.again[4], 3), rep(denom.again[5], 3),
          rep(denom.again[6], 3) )

tEmp<- rep(0, 18)
for(a in 1:18){
  tEmp[a]<- p.ytilde[a, 1]*lit.p[a,2]/temP[a]
}

tEmp<- data.frame(final.p=tEmp, age = p.ytilde$ages, diagnosis=p.ytilde$diagnosis)
tEmp
# ---- plot

tEmp[is.na(tEmp)]<- 2
tEmp

####################################################
# Plots similar to Figure 5 for each quantile 

windows()
ggplot(data=tEmp, aes(x = age, y=final.p, colour=diagnosis)) + geom_point() + geom_line() +
  ggtitle(paste("SIMULATED ", quant)) + ylab("Probability") + xlab("Age") +theme_bw() + theme(legend.position = "none")
# Contact primary author for details on Monte Carlo error intervals
#------------------------------------------------------------------------------------------------------------------------

# **************************  third quantile 25.75
top.q<- hippo.ranges[c(3,7,11,15,19,23),]
rownames(top.q)<- seq(1:6)
top.q
quant<- "25 - 50th Quantile"

#----------
# Compute from the model: P(y_tilde| diagnosis, age)
temp<- data.frame(matrix(data=rep(0, no.chains*6*3), nrow = no.chains*3, ncol=6))

for(j in 1:6){  
  for(u in 1:no.chains*3){   
    if(top.q[j,2] < Age.chains2[u,j] & Age.chains2[u,j] < top.q[j,3]){
      temp[u,j]<- 1} 
  }
}
#head(temp)
#colSums(temp)
temp$Diag<- c(rep("HC", no.chains), rep("MCI", no.chains), rep("AD", no.chains))

# now see what each diagnosis comes as
hc.p<- subset(temp, Diag == "HC")
all.ages.hc.p<- colSums(hc.p[, 1:6])/no.chains

mci.p<- subset(temp, Diag == "MCI")
all.ages.mci.p<- colSums(mci.p[1:6])/no.chains

ad.p<- subset(temp, Diag == "AD")
all.ages.ad.p<- colSums(ad.p[1:6])/no.chains

# put it all into one data frame
p.ytilde.Diag.age<- data.frame(probs.mdl = c(all.ages.hc.p, all.ages.mci.p, all.ages.ad.p), ages = rep( c(60,65,70,75,80,85), 3),
                               diagnosis = c(rep("HC", 6), rep("MCI", 6), rep("AD", 6) ) )   
p.ytilde.Diag.age
p.ytilde<- p.ytilde.Diag.age[order(p.ytilde.Diag.age$ages),]
p.ytilde

#---- now to include probabilities from literature, expression for denominator

denom.again<- c(0,0,0,0,0,0)
d<- c(3,6,9,12,15,18)
for(i in 1:6){
  teMp<- d[i]
  denom.again[i] <-  p.ytilde[teMp-2,1]*lit.p[teMp-2,2] +  p.ytilde[teMp-1,1]*lit.p[teMp-1,2] + p.ytilde[teMp,1]*lit.p[teMp,2] 
}
denom.again 

#----- now estimate final probs

temP<- c( rep(denom.again[1], 3), rep(denom.again[2], 3), rep(denom.again[3], 3), rep(denom.again[4], 3), rep(denom.again[5], 3),
          rep(denom.again[6], 3) )

tEmp<- rep(0, 18)
for(a in 1:18){
  tEmp[a]<- p.ytilde[a, 1]*lit.p[a,2]/temP[a]
}

tEmp<- data.frame(final.p=tEmp, age = p.ytilde$ages, diagnosis=p.ytilde$diagnosis)
tEmp
# ---- plot

tEmp[is.na(tEmp)]<- 2
tEmp

####################################################
# Plots similar to Figure 5 for each quantile 

windows()
ggplot(data=tEmp, aes(x = age, y=final.p, colour=diagnosis)) + geom_point() + geom_line() +
  ggtitle(paste("SIMULATED ", quant)) + ylab("Probability") + xlab("Age") +theme_bw() + theme(legend.position = "none")
# Contact primary author for details on Monte Carlo error intervals
#------------------------------------------------------------------------------------------------------------------------

# *****************fourth quantile
top.q<- hippo.ranges[c(4,8,12,16,20,24),]
rownames(top.q)<- seq(1:6)
top.q
quant<- "15 - 25th Quantile"

#----------
# Compute from the model: P(y_tilde| diagnosis, age)
temp<- data.frame(matrix(data=rep(0, no.chains*6*3), nrow = no.chains*3, ncol=6))

for(j in 1:6){  
  for(u in 1:no.chains*3){   
    if(top.q[j,2] < Age.chains2[u,j] & Age.chains2[u,j] < top.q[j,3]){
      temp[u,j]<- 1} 
  }
}
#head(temp)
#colSums(temp)
temp$Diag<- c(rep("HC", no.chains), rep("MCI", no.chains), rep("AD", no.chains))

# now see what each diagnosis comes as
hc.p<- subset(temp, Diag == "HC")
all.ages.hc.p<- colSums(hc.p[, 1:6])/no.chains

mci.p<- subset(temp, Diag == "MCI")
all.ages.mci.p<- colSums(mci.p[1:6])/no.chains

ad.p<- subset(temp, Diag == "AD")
all.ages.ad.p<- colSums(ad.p[1:6])/no.chains

# put it all into one data frame
p.ytilde.Diag.age<- data.frame(probs.mdl = c(all.ages.hc.p, all.ages.mci.p, all.ages.ad.p), ages = rep( c(60,65,70,75,80,85), 3),
                               diagnosis = c(rep("HC", 6), rep("MCI", 6), rep("AD", 6) ) )   
p.ytilde.Diag.age
p.ytilde<- p.ytilde.Diag.age[order(p.ytilde.Diag.age$ages),]
p.ytilde

#---- now to include probabilities from literature, expression for denominator

denom.again<- c(0,0,0,0,0,0)
d<- c(3,6,9,12,15,18)
for(i in 1:6){
  teMp<- d[i]
  denom.again[i] <-  p.ytilde[teMp-2,1]*lit.p[teMp-2,2] +  p.ytilde[teMp-1,1]*lit.p[teMp-1,2] + p.ytilde[teMp,1]*lit.p[teMp,2] 
}
denom.again 

#----- now estimate final probs

temP<- c( rep(denom.again[1], 3), rep(denom.again[2], 3), rep(denom.again[3], 3), rep(denom.again[4], 3), rep(denom.again[5], 3),
          rep(denom.again[6], 3) )

tEmp<- rep(0, 18)
for(a in 1:18){
  tEmp[a]<- p.ytilde[a, 1]*lit.p[a,2]/temP[a]
}

tEmp<- data.frame(final.p=tEmp, age = p.ytilde$ages, diagnosis=p.ytilde$diagnosis)
tEmp
# ---- plot

tEmp[is.na(tEmp)]<- 2
tEmp

####################################################
# Plots similar to Figure 5 for each quantile 

windows()
ggplot(data=tEmp, aes(x = age, y=final.p, colour=diagnosis)) + geom_point() + geom_line() +
  ggtitle(paste("SIMULATED ", quant)) + ylab("Probability") + xlab("Age") +theme_bw() + theme(legend.position = "none")
# Contact primary author for details on Monte Carlo error intervals
#------------------------------------------------------------------------------------------------------------------------
