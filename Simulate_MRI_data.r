#############################################
# Simulate Ventricle.ICV and Sum of left and right Hippocamputs ICV 
# MRI values: to replicate results in manuscript
#
# Comparisons of neurodegeneration over time between healthy ageing and Alzheimer's 
# disease cohorts
#
# *** insert ePubs or journal reference ***

# See supplementary material for posterior mean parameter estimates
# for model (3) in manuscript.

library(ggplot2)
  
N<- 260  # 260 participants
tp<- 4 # number of maximum timepoints observed for each person
# note: StndAge must be sequential. I.e -2, -1.5, 0, 1.
#                                  not  -2, 1, -1.5, 0.

# StndAge wide
baseline.age<- runif(N, min = -2, max = 1.2)
age.t2<- age.t3<- age.t4<- c()
for(j in 1:N){
  age.t2[j]<- baseline.age[j] + 0.2  # StndAge of 0.2 is approximately 1.5 years
  age.t3[j]<- age.t2[j] + 0.2
  age.t4[j]<- age.t3[j] + 0.2
}
age<- data.frame(X1 = baseline.age, X2=age.t2, X3 = age.t3, X4 = age.t4)

# now create missing values -> to simulate unbalanced data
# recall in longitudinal neuroimaging studies patient drop out as 
# time progressess increases

seed.v<- 2 # change to generate different data 
set.seed(seed.v)
tp1.miss<- sample(1:N, size= 20)
set.seed(seed.v+1)
tp2.miss<- sample(1:N, size = 50)
set.seed(seed.v+2)
tp3.miss<- sample(1:N, size = 70)
set.seed(seed.v+3)
tp4.miss<- sample(1:N, size = 80)

age[tp1.miss, 1] <- NA # first, second, third and fourth columns
age[tp2.miss, 2]<- NA  # insert missing values
age[tp3.miss, 3]<- NA
age[tp4.miss, 4]<- NA

# create diagnosis wide data - with 27 converters.
# Figure 4 of manuscript describes cohort with 157 HC, 34 MCI and 42 AD participants

# 27 participants converted throughout the study: 
sim.diag<- c(rep("HC", 157*tp), rep("MCI", 34*tp), rep("AD", 42*tp), 
             rep(c("HC","HC", "MCI", "MCI"), times = 8), # 8 converted HC -> MCI
             rep(c("HC", "AD", "AD", "AD"), times = 4), # say 4 converted HC -> AD
             rep(c("MCI", "MCI", "AD", "AD"), times = 15)) # and say 15 converted from MCI -> AD

diag<- data.frame(matrix(0, nrow = N, ncol = tp) )

count<- 1
for(r in 1:N){
  for(c in 1:tp){
    diag[r,c]<- sim.diag[count]
    count<- count+1
  }
}

diag[tp1.miss, 1]<-NA
diag[tp2.miss, 2]<-NA
diag[tp3.miss, 3]<-NA
diag[tp4.miss, 4]<-NA

diag  # check that everyone has AT LEAST one observation

#save(diag, file = "diag.Rdata")
final.sim.dat<- data.frame(ID = rep(1:260, times = 4), Diagnosis=melt(diag, measure.vars=c("X1", "X2", "X3", "X4")), 
                           StndAge = melt(age, measure.vars= c("X1", "X2", "X3", "X4")) )
head(final.sim.dat)

final.sim.dat<- final.sim.dat[, -c(2, 4)]
colnames(final.sim.dat)<- c("ID", "Diagnosis", "StndAge")
final.sim.dat<- na.omit(final.sim.dat)  # NA removed and long format

final.sim.dat<- final.sim.dat[order(final.sim.dat$ID),]
rownames(final.sim.dat)<- seq(1:dim(final.sim.dat)[1])

final.sim.dat$Diagnosis<- factor(final.sim.dat$Diagnosis, levels = c("HC", "MCI", "AD"))
contrasts(final.sim.dat$Diagnosis) # Unfortunately in BUGS and rJAGS we need to define our own contrasts/
                                   # binarise the categorical variable 
final.sim.dat$mci.bin <-ifelse(final.sim.dat$Diagnosis == "MCI", 1, 0)
final.sim.dat$ad.bin <-ifelse(final.sim.dat$Diagnosis == "AD", 1, 0)

   #####################################################################

# create random effects for HC, MCI and AD participants

                # *************
                # Ventricle   *
                # *************

# from model (3) in the manuscript, posterior means for the parameters are

# ventricle ICV
b0.v<- 2.364012e-02; b1.v<- 5.689918e-03
b4.v<- 9.371770e-04; b5.v<- 3.864530e-03

b2.v<- 9.208147e-04; b3.v<- 5.102918e-03; s2.v<- 2.298996e-06 
S.v<- matrix(c(1.384803e-04, 2.411531e-05, -2.591587e-06, -1.308706e-05,
               2.411531e-05, 9.557799e-06, -5.639764e-07, -1.617011e-06,
              -2.591587e-06, -5.639764e-07, 1.019274e-05, 1.229625e-06,
              -1.308706e-05, -1.617011e-06, 1.229625e-06, 3.339590e-05), nrow=4, ncol=4)

# Sum of left and right hippo ICV
b0.h<- 4.024719e-01; b1.h<- -1.293484e-02 
b4.h<- -8.614252e-03; b5.h<- -1.482345e-02 

b2.h<- -1.096659e-02; b3.h<- -2.288477e-02; s2.h<- 4.228257e-05
S.h<- matrix(c(1.166047e-03, 2.857232e-04, -1.602743e-04, -3.006436e-04,
               2.857232e-04, 1.223680e-04, 3.145111e-05, 2.769425e-06,
              -1.602743e-04, 3.145111e-05, 3.271362e-04, 3.211273e-04,
              -3.006436e-04, 2.769425e-06, 3.211273e-04, 4.777180e-04), nrow=4, ncol=4)

Vent.ICV<- c(); beta.v<- beta.h<- matrix(0, nrow = 260, ncol = 4)
Hippo.ICV<- c(); 
mu.h<- c(b0.h, b1.h, b4.h, b5.h)
mu.v<- c(b0.v, b1.v, b4.v, b5.v)
person<- as.numeric(as.character(final.sim.dat$ID))

for(n in 1:N){  # random effects
  beta.v[n, 1:4] <- mvrnorm(1, mu=mu.v, Sigma=S.v)
  beta.h[n, 1:4]<- mvrnorm(1, mu=mu.h, Sigma=S.h)
}

# ********************************************
for(i in 1:dim(final.sim.dat)[1]){
  age<- final.sim.dat$StndAge[i]
  x1<- final.sim.dat$mci.bin[i]
  x2<- final.sim.dat$ad.bin[i]

# generate Ventricle values  
  y.v<- beta.v[person[i], 1] + b2.v*x1 +b3.v*x2 +(beta.v[person[i], 2] + beta.v[person[i], 3]*x1 +
                                                         beta.v[person[i], 4]*x2)*age
  Vent.ICV[i]<- rnorm(1, y.v, sqrt(s2.v))

# generate Sum left and right hippocampus values
  y.h<- beta.h[person[i], 1] + b2.h*x1 +b3.h*x2 +(beta.h[person[i], 2] + beta.h[person[i], 3]*x1 +
                                                         beta.h[person[i], 4]*x2)*age
  Hippo.ICV[i]<- rnorm(1, y.h, sqrt(s2.h))
}

final.sim.dat$Vent.ICV<- Vent.ICV
final.sim.dat$Hippo.ICV<- Hippo.ICV

save(final.sim.dat, file = "final.sim.dat.Rdata")

                  # ****************************************** 
                  # visualise ventricle spaghetti plot data  *
                  # ******************************************

final.sim.dat$ID<- factor(final.sim.dat$ID)
final.sim.dat$Diagnosis<- factor(final.sim.dat$Diagnosis)

windows() # VENTRICLE
ggplot(final.sim.dat, aes(x = StndAge, y = Vent.ICV, group = ID, colour = Diagnosis)) + geom_line()  + 
  ggtitle("Simulated Ventricle.ICV MRI data") + theme_bw() #+theme(legend.position="none") 
# note there are negative values here, and this is not just because of the s2.v

# mebe remove those people who has at least one negative MRI value
#ppl.neg.vol<- unique(subset(final.sim.dat, Vent.ICV < 0)$ID)
#f.sim.d2<- subset(final.sim.dat, !ID %in%ppl.neg.vol)

#f.sim.d2$ID<- factor(as.character(f.sim.d2$ID))
#windows()
#ggplot(f.sim.d2, aes(x = StndAge, y = Vent.ICV, group = ID, colour = Diagnosis)) + geom_line()  + 
#  ggtitle("Corrrected Simulated Ventricle.ICV MRI.ICV data") + theme_bw()

                    # ******************************
windows() # HIPPO
ggplot(final.sim.dat, aes(x=StndAge, y=Hippo.ICV, group =ID, colour = Diagnosis)) + geom_line() +
  ggtitle("Simulated sum of left and right hippocampus MRI.ICV data\n scaled by 100") + theme_bw()
