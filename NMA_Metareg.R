
#packages
require(meta)
require(netmeta)
require(readxl)
require(dplyr)
require(tidyr)
require(magrittr)
require(rjags)
require(gemtc)
require(dmetar)

#load data
dat<-CollatedData_2020_10_09_same_int_removed_imputed_MBNMA <- read_excel("CollatedData_2020_10_09_same int_removed_imputed_MBNMA.xlsx")



#select columns for pain short-term
dat1 <-
  dat %>%
  select(INFO_Author.Year, COMPARATOR_SHORT_HAND,Back.pain.or.Pain_Is.more.or.less.better.,
         Back.Pain...1d.but..3mo_Mean, Back.Pain...1d.but..3mo_SD, Back.Pain...1d.but..3mo_N,Baseline.Back.Pain_Mean,Baseline.Back.Pain_SD,POPULATION_N.enrolled)
dat1


#drop NAs
dat1 <- dat1 %>% drop_na()
dat1



#Calculate SMDs
p1<-pairwise(treat=COMPARATOR_SHORT_HAND,
             mean=Back.Pain...1d.but..3mo_Mean,sd=Back.Pain...1d.but..3mo_SD,n=Back.Pain...1d.but..3mo_N,
             data =dat1, studlab=INFO_Author.Year,sm="SMD")
p1

#rescale -/+1
p1 <- transform(p1, p1$TE==ifelse(Back.pain.or.Pain_Is.more.or.less.better.=="more.is.better", p1$TE*-1, p1$TE*1))
p1


#calculate Baseline normalized mean

#calculate standard errors of each group
se1 <- p1$Baseline.Back.Pain_SD1/(p1$POPULATION_N.enrolled^0.5)
se2 <- p1$Baseline.Back.Pain_SD2/(p1$POPULATION_N.enrolled^0.5)

#calculate inverse variance weights for each group
w1 <- 1/(se1^2)
w2 <- 1/(se2^2)


#pooled mean
M<-(p1$Baseline.Back.Pain_Mean1*w1+p1$Baseline.Back.Pain_Mean2*w2)/(w1 + w2)


#pooled variance of the mean

V<-1/(w1+w2)

#pooled SE
SE <- V^0.5

#pooled SD
SD<-SE * ((p1$POPULATION_N.enrolled*2)^0.5)

#normalized mean
NM<-M/SD

p1<-data.frame(p1,NM)



#select columns for pain (short-term) effect sizes
p1 <-

  p1 %>% select(studlab, treat1, treat2, TE, seTE, NM)
p1


#transform for getmc data format
p1 %>%
  dplyr::select(1:6) %>%
  pivot_longer(-studlab,
               names_to = c(".value"),
               names_pattern = "(..)") %>%
  set_colnames(c("study", "treatment",
                 "diff", "std.err","NM")) -> data


#select columns for pain (short-term) effect sizes
data1 <-
  data %>%
  select(study,treatment,diff,std.err)
data1

data1<-as.data.frame(data1)


#################################BAYESIAN NMA###########################################################
#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/bayesnma.html

#network graph
network <- mtc.network(data.re = data1)
summary(network)
plot(network)


#remove rows ---no closed network
data1 <- data1[-c(7,8),]

# Just to check here - dataset is given as SMD treatment contrasts right?
#network
network <- mtc.network(data.re = data1)
plot(network)

#build model and run simulations

model <- mtc.model(network,
                   linearModel = "random",
                   n.chain = 4)
mcmc1 <- mtc.run(model, n.adapt = 50, n.iter = 1000, thin = 10) ####take #2
mcmc2 <- mtc.run(model, n.adapt = 5000, n.iter = 100000, thin = 10)

#Model diagnostics
gelman.plot(mcmc1)
gelman.plot(mcmc2)

gelman.diag(mcmc1)$mpsrf
gelman.diag(mcmc2)$mpsrf

# Based on the bgrplots (gelman.plot) I would use more adapt samples - can probs use fewer monitored iterations too
mcmc3 <- mtc.run(model, n.adapt = 10000, n.iter = 20000, thin = 10)
gelman.plot(mcmc3)
gelman.diag(mcmc3)$mpsrf

#Assessing inconsistency
nodesplit <- mtc.nodesplit(network,
                           linearModel = "random",
                           n.adapt = 10000,
                           n.iter = 20000,
                           thin = 10)

# Some really quite substantial inconsistency here which would be cause for concern
# Particularly for PT vs Exercise and PT vs Pharma where the direct and indirect suggest opposite effects
# Probably due to lumping of treatments - assuming all pharma treatments to be equivalent seems like a very strong assumption to me
# The same perhaps goes for other treatments in the network too, but I don't have clinical expertise to be able to know about this
summary(nodesplit)
plot(summary(nodesplit))

#Results
rank <- rank.probability(mcmc2, preferredDirection = -1)
plot(rank, beside=TRUE, cex.names=0.5)

rank.probability <- rank.probability(mcmc2)
sucra(rank.probability, lower.is.better = TRUE)
results <- relative.effect.table(mcmc2)

#####################################META-REGRESSIION#######################################################

#create new data frame
data<-as.data.frame(data)
#remove rows ---no closed network
data <- data[-c(7,8),]

covariate<-data.frame(data$study,data$NM)

#rename columns
names(covariate) <- c("study", "NM")


network.mr <- mtc.network(data.re = data1,
                          studies = covariate)

regressor <- list(coefficient = 'shared',    #####I am not 100% clear about the specification
                  variable = 'NM',
                  control = 'Placebo')

model.mr <- mtc.model(network.mr,
                      type = "regression",
                      regressor = regressor)

mcmc3 <- mtc.run(model.mr,
                 n.adapt = 10000,
                 n.iter = 20000,
                 thin = 10)
summary(mcmc3)
####The results for our covariate are reported under B####
####median=-0.150369 95%CrI (-1.9580, 1.7036) includes zero


# So above you assume a shared regression effect (that the regressor has the same on all treatment effects)
# We could relax this assumption by assuming an exchangeable regressor effect
# The difference between these is akin to a fixed vs random effects meta-analysis, except that the fixed/random effects are
#on the regressor effect and are labelled as shared/exchangeable.
regressor <- list(coefficient = 'exchangeable',    #####I am not 100% clear about the specification
                  variable = 'NM',
                  control = 'Placebo')

model.mr <- mtc.model(network.mr,
                      type = "regression",
                      regressor = regressor)

mcmc4 <- mtc.run(model.mr,
                 n.adapt = 10000,
                 n.iter = 20000,
                 thin = 10)
summary(mcmc4)

# The most complex model with fewest assumptions would be to assume that the regressor effects are "unrelated" between
#different treatments. But there will probably not be enough information to assume this - you would need multiple studies
#of each treatment comparison. A "middle ground" approach would be to assume a shared effect within different classes, which
#allows for a bit more flexibility. So for treatments investigated in fewer studies you could assume it had the same class
#as for another treamtent.
# This is probably something that we need to sort ASAP actually - which classes are we going to assume the different treatments
#fall into?


###Use DIC (deviance information criterion) to compare models.Lower DIC values mean better fit.

summary(mcmc4)$DIC
summary(mcmc3)$DIC
summary(mcmc2)$DIC
###here almost the same results.Imho there is no influence of Baseline.pain in the short term.
# You are correct in your conclusion for model 3 vs model 2....but if you look at mcmc4 there is a big reduction in DIC
# Also there is a reduction in the between-study SD, suggesting that including this regressor explains some of the heterogeneity
# This suggests that for some treatments it is an important regressor, whilst for others it may not be

# If we look at the model results:
summary(mcmc4)
# beta[10] has a very negative coefficient
# beta[2] does not have such a negative coefficient

# So this suggests that for treatment 10 having a lower baseline score substantially reduces the treatment effect, whereas this
#reduction is much less for treatment 2
# Note that these coefficients will be "shunk" towards the mean exchangeable regressor effect (B)

