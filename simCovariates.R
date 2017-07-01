# Structural topic models -------------------------------------------------

# Structural topic models allow us to model topic prevalence with a set of 
# covariates of interest. In what follows, a data generating function similar to
# the approach displayed in the main document is used to generate the data 
# assuming the effects of two covariates. See the simTreatmentEffect.R for the a
# model with the group effect alone.

# Packages needed: DirichletReg, topicmodels, stm, heatmaply, tidyverse, mgcv

# See http://www.structuraltopicmodel.com/ for more info and references.


# Data generating function and data creation ------------------------------

gendat = function(nt=3, nd=500, nw=40){
  # nt: number of topics
  # nd: number of docs
  # nw: average number of words per doc
  
  # group: a binary group factor covariate
  # x: a continuous covariate
  # B: topic is a distribution over the fixed vocabulary; e.g. apple = .02, berry = .01, cat = .0001
  # theta: topic proportions for each document where theta[i,] is the topic distribution for the ith document
  # Z: term probabilities given topics
  # wordList: are the observed words
  #   
  # require(DirichletReg)  # for vectorized rdirichlet
  # x = sort(rnorm(nd))    # if you want a more stark visual for topic probs
  x = rnorm(nd)
  group = rep(0:1, nd/2)
  X = cbind(1, group, x)
  t1 = X %*% c(0,  1,  2)  # positive effects of group and treatment for topic 1
  t2 = X %*% c(0, -1, -2)  # negative effects of group and treatment for topic 2
  t3 = 0                   # no effects for topic 3
  alpha = plogis(cbind(t1, t2, t3))
  
  Nd = rpois(nd, 40)
  G0 = 1/3
  theta = DirichletReg::rdirichlet(nd, alpha * G0)  # lda approach
  B = DirichletReg::rdirichlet(n=nw, alpha=rep(.05, nt)) 
  Z = tcrossprod(theta, B)
  wdList = vector('list', nd)
  for (i in 1:nrow(Z))  wdList[[i]] = t(rmultinom(1, Nd[i], Z[i,]))
  ldaform = sapply(wdList, function(x) rbind(1:40,x), simplify = F)
  wd = do.call(rbind, wdList)
  wdList = lapply(wdList, function(wds) rep(paste0('word',1:length(wds)), wds))
  require(stringr)
  
  return(list(dtmat=wd, BagofWords=wdList, ldaDocFormat=ldaform, X=as.data.frame(X[,-1])))
}


# debugonce(gendat)
# gendat()

set.seed(4377)
ndocs = 1000
docs = gendat(nd=ndocs)


# Standard LDA ------------------------------------------------------------

library(topicmodels)
lda_0 = posterior(LDA(docs[[1]], k=3, method = 'VEM'))  

#§ Explore LDA with no covariates ----

library(heatmaply)
heatmaply(lda_0$topics, Rowv = NA, Colv = NA, colors=viridis::inferno(100),
          showticklabels=FALSE, plot_method='plotly', fontsize_row=0, fontsize_col=0)


library(tidyverse)
lm_prob  = lm(lda_0$topics ~., data=docs$X)  # coefs on probability scale, will be similar to stm
lm_logit = lm(qlogis(lda_0$topics) ~., data=docs$X)  # more reflect the coefs used in the function
lm_prob
lm_logit

# the following models would treat the topics more on their terms via beta or dirichlet regression.
library(mgcv); library(DirichletReg)
beta_reg = apply(lda_0$topics, 2, function(y) gam(y ~ group + x, data=docs$X, family='betar'))
sapply(beta_reg, coef)
df = docs$X
df$topics = DR_data(lda_0$topics)
dir_reg = DirichReg(topics ~ group + x, data=df, control=list(tol1=1e-12, tol2=1e-12))
coef(dir_reg)


gdat_orig = data.frame(docs$X, lda_0$topics) %>% 
  gather(key=topic, value=proportion, X1, X2, X3)
gdat_fits = data.frame(docs$X, X=fitted(dir_reg)) %>% 
  dplyr::rename(X1 = X.1, X2 = X.2, X3 = X.3) %>% 
  gather(key=topic, value=fitted, X1, X2, X3)
gdat = left_join(gdat_orig, gdat_fits) %>% 
  mutate(group=factor(group))
head(gdat)

gdat %>%
  group_by(topic, group) %>%
  arrange(desc(x)) %>%
  plot_ly(color=~topic) %>%
  add_lines(x=~x, y=~fitted, linetype=~group, colors=viridis::viridis(3)) %>% 
  lazerhawk::theme_plotly()


LDA_topic_props = data.frame(group=docs$X$group, lda_0$topics) %>% 
  group_by(group) %>% 
  summarise_all(mean)
LDA_topic_props

LDA_diffs = LDA_topic_props %>% 
  select(-group) %>%
  sapply(diff) %>%
  sort %>%
  round(3)
LDA_diffs


# Structured topic models with stm package --------------------------------

library(stm)
stm_1 = stm(docs[[3]], vocab=paste0('word',1:40), K=3, prevalence = ~ group + x, data=docs$X)

heatmaply(stm_1$theta, Rowv = NA, Colv = NA, colors=viridis::inferno(100),
          showticklabels=FALSE, plot_method='plotly',
          fontsize_row=0, fontsize_col=0)


#§ Use stm to estimate the effect directly and visualize ----

# match topics; labels are arbitrary, so a cor of ~1 would indicate the same topic
cor(data.frame(lda=lda_0$topics), data.frame(stm=stm_1$theta))

# estimate effects
ee = estimateEffect(~group+x, stm_1, docs$X)
summary(ee)
plot(ee, 'group',  method = 'difference', cov.value1='treatment', cov.value2 = 'control')
plot(ee, 'x',  method = 'continuous')
