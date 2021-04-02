
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

#Assessing inconsistency
nodesplit <- mtc.nodesplit(network,
                           linearModel = "random",
                           n.adapt = 5000,
                           n.iter = 100000,
                           thin = 10)

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
                 n.adapt = 5000,
                 n.iter = 100000,
                 thin = 10)
summary(mcmc3)
####The results for our covariate are reported under B####
####median=-0.150369 95%CrI (-1.9580, 1.7036) includes zero


###Use DIC (deviance information criterion) to compare models.Lower DIC values mean better fit.

summary(mcmc3)$DIC
summary(mcmc2)$DIC
###here almost the same results.Imho there is no influence of Baseline.pain in the short term.
