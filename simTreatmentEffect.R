
# Structural topic models -------------------------------------------------

# Structural topic models allow us to model topic prevalence with a set of
# covariates of interest. In what follows, a data generating function similar to
# the approach displayed in the main document is used to generate the data
# assuming a binary grouping of documents (e.g. perhaps they come from two
# different sources). The default makeup suggests a negative group difference
# for the first topic ('control' > 'treatment' group), no difference for the
# second, and a positive difference (same magnitude) for the third topic. We
# then will use standard LDA and STM to explore the group differences.

# Packages needed: DirichletReg, topicmodels, stm, heatmaply, tidyverse

# See http://www.structuraltopicmodel.com/ for more info and references.



# Data generating function and data creation ------------------------------

gendat = function(nt=3, nd=500, nw=40){
  # nt: number of topics
  # nd: number of docs
  # nw: average number of words per doc

  # B: topic is a distribution over the fixed vocabulary; e.g. apple = .02, berry = .01, cat = .0001
  # theta: topic proportions for each document where theta[i,] is the topic distribution for the ith document
  # Z: term probabilities given topics
  # wordList: are the observed words

  # Preliminary setup
  # require(DirichletReg)  # for vectorized rdirichlet
  alpha0 = c(.3,.4,.3)     # topic props control
  alpha1 = c(.1,.4,.5)     # topic props treatment group +-.2 for two of the topics
  treatment = rbind(t(replicate(nd/2, alpha0)),
                    t(replicate(nd/2, alpha1)))
  Nd = rpois(nd, 40)
  G0 = 1/3
  
  # generate documents
  theta = DirichletReg::rdirichlet(nd, treatment * G0)  # note that topic prevalence depends on group
  B = DirichletReg::rdirichlet(n=nw, alpha=rep(.05, 3)) 
  Z = tcrossprod(theta, B)
  wordList = vector('list', nd)
  for (i in 1:nrow(Z))  wordList[[i]] = t(rmultinom(1, Nd[i], Z[i,]))
  ldaform = sapply(wordList, function(x) rbind(1:40, x), simplify = F) # for use with stm
  
  # cleanup
  wd = do.call(rbind, wordList)
  wordList = lapply(wordList, function(wds) rep(paste0('word', 1:length(wds)), wds))
  
  return(list(dtmat=wd, BagofWords=wordList, ldaDocFormat=ldaform))
}

ndocs= 500
treatment = factor(rep(0:1, e=ndocs/2), labels = c('control', 'treatment'))
docs = gendat(nd=ndocs)
str(docs, 1)

# Standard LDA ------------------------------------------------------------

library(topicmodels)
lda_0 = posterior(LDA(docs[[1]], k=3, method = 'VEM', control=list(nstart=10, verbose=100)))  


#§ Explore ----

library(heatmaply)
heatmaply(lda_0$topics, Rowv = NA, Colv = NA, colors=viridis::magma(100),
          showticklabels=FALSE, plot_method='plotly', fontsize_row=0, fontsize_col=0)

library(tidyverse)
LDA_topic_props = data.frame(treatment=treatment, lda_0$topics) %>% 
  group_by(treatment) %>% 
  summarise_all(mean)
LDA_topic_props

LDA_diffs = LDA_topic_props %>% 
  select(-treatment) %>%
  sapply(diff) %>%
  sort %>%
  round(3)
LDA_diffs



# Structured topic models with stm package --------------------------------

library(stm)
stm_0 = stm(docs[[3]], vocab=paste0('word', 1:40), K=3) # this is just LDA
stm_1 = stm(docs[[3]], vocab=paste0('word', 1:40), K=3, 
            prevalence = ~treatment,
            data=data.frame(treatment))

heatmaply(stm_0$theta, Rowv = NA, Colv = NA, colors=viridis::magma(100),
          showticklabels=FALSE, plot_method='plotly', 
          fontsize_row=0, fontsize_col=0)

#§ Explore STM with no covariate (same as lda_0) ----

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

#§ Explore STM with covariate ----

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



#§ Use stm to estimate the effect directly and visualize ----

ee = estimateEffect(~treatment, stm_1, data.frame(treatment))
summary(ee)
plot(
  ee,
  covariate  = 'treatment',
  method     = 'difference',
  cov.value1 = 'treatment',
  cov.value2 = 'control'
)



# Compare all results -----------------------------------------------------

comparison = rbind(LDA_diffs, stm_0_diffs, stm_1_diffs)
colnames(comparison) = c('Topic 1', 'Topic 2', 'Topic 3') # reset arbitrary topic names
comparison

