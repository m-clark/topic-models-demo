### Stm appendix to open-ended survey responses

gendat = function(nt=3, nd=500, nw=40){
#   B_k topic is a distribution over the fixed vocabulary; e.g. apple = .02, berry = .01, cat = .0001
#   theta_d topic proportions for dth document where theta_dk is the topic prop for the kth topic in doc d
#   z_d topic assignment for dth doc where z_dn is topic ass for nth word in doc d
#   w_d are the observed words in doc d where w_dn is the nth word in doc d
#   
  # require(DirichletReg) # for vectorized rdirichlet
  B_k = DirichletReg::rdirichlet(n=nw, alpha=rep(.05,3)) 
  alpha0 = c(.3,.4,.3)            # topic props control
  alpha1 = c(.1,.4,.5)            # topic props treatment group +-.2 for two of the topics
  treatment = rbind(t(replicate(nd/2, alpha0)),
                    t(replicate(nd/2, alpha1)))
  Nd = rpois(nd, 40)
  G0 = 1/3
  theta_dk = DirichletReg::rdirichlet(nd, treatment * G0)
  thetabeta = tcrossprod(theta_dk, B_k)
  wdList = vector('list', nd)
  for (i in 1:nrow(thetabeta))  wdList[[i]] = t(rmultinom(1, Nd[i], thetabeta[i,]))
  ldaform = sapply(wdList, function(x) rbind(1:40,x), simplify = F)
  wd = do.call(rbind, wdList)
  wdList = lapply(wdList, function(wds) rep(paste0('word',1:length(wds)), wds))
  require(stringr)
  
  return(list(dtmat=wd, BagofWords=wdList, ldaDocFormat=ldaform))
}
# 
# debugonce(gendat)

ndocs= 500
treatment = factor(rep(0:1, e=ndocs/2), labels = c('control', 'treatment'))
docs = gendat(nd=ndocs)

# # for reordering
# idx = sample(1:ndocs)
# treatment = treatment[idx]
#   
# str(docs, 1)
# str(docs[[2]][[1]])
# (docs[[3]][[1]])


### LDA
library("topicmodels")
testLDA = posterior(LDA(docs[[1]], k=3, method = 'VEM', control=list(nstart=10, verbose=100)))  
# testLDA = posterior(LDA(docs[[1]][idx,], k=3, method = 'VEM'))  

library(heatmaply)
heatmaply(testLDA$topics, Rowv = NA, Colv = NA, colors=viridis::magma(100),
          showticklabels=FALSE, plot_method='plotly', fontsize_row=0, fontsize_col=0)

library(tidyverse)
LDA_topic_props = data.frame(treatment=treatment, testLDA$topics) %>% 
  group_by(treatment) %>% 
  summarise_all(mean)
LDA_topic_props

LDA_diffs = LDA_topic_props %>% 
  select(-treatment) %>%
  sapply(diff) %>%
  sort %>%
  round(3)
LDA_diffs


### STM

library(stm)
stm_0 = stm(docs[[3]], vocab=paste0('word', 1:40), K=3)
stm_1 = stm(docs[[3]], vocab=paste0('word', 1:40), K=3, prevalence = ~treatment, 
            data=data.frame(treatment))

heatmaply(stm_0$theta, Rowv = NA, Colv = NA, colors=viridis::magma(100),
          showticklabels=FALSE, plot_method='plotly', 
          fontsize_row=0, fontsize_col=0)

stm_0_topic_props = data.frame(treatment=treatment, stm_0$theta) %>% 
  group_by(treatment) %>% 
  summarise_all(mean)

stm_0_topic_props

stm_0_diffs = stm_0_topic_props %>% 
  select(-treatment) %>%
  sapply(diff) %>%
  sort %>%
  round(3)

stm_0_diffs

# STM with covariate

heatmaply(stm_1$theta, Rowv = NA, Colv = NA, colors=viridis::magma(100),
          showticklabels=FALSE, plot_method='plotly', 
          fontsize_row=0, fontsize_col=0)

stm_1_topic_props = data.frame(treatment=treatment, stm_1$theta) %>% 
  group_by(treatment) %>% 
  summarise_all(mean)

stm_1_topic_props

stm_1_diffs = stm_1_topic_props %>% 
  select(-treatment) %>%
  sapply(diff) %>%
  sort %>%
  round(3)

stm_1_diffs




ee = estimateEffect(~treatment, stm_1, data.frame(treatment))
plot(ee, covariate='treatment',  method = 'difference', cov.value1='control', cov.value2 = 'treatment')

rbind(LDA_diffs, stm_0_diffs, stm_1_diffs)


