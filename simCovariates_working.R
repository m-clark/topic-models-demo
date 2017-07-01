# Structural topic models -------------------------------------------------

# Structural topic models allow us to model topic prevalence with a set of 
# covariates of interest. In what follows, a data generating function similar to
# the approach displayed in the main document is used to generate the data 
# assuming the effects of two covariates. See the simTreatmentEffect.R for the a
# model with the group effect alone.

# See http://www.structuraltopicmodel.com/ for more info and references.



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
  # require(DirichletReg)
  x = sort(rnorm(nd))
  group = rep(0:1, nd/2)
  X = cbind(1, group, x)
  t1 = X %*% c(0, 1, 2)
  t2 = X %*% c(0, -1, -2)
  t3 = .4
  alpha = plogis(cbind(t1, t2, t3))
  
  Nd = rpois(nd, 40)
  G0 = 1/3
  theta = rdirichlet(nd, alpha * G0)  # lda approach
  # theta = plogis(cbind(t1,t2,0))    # logistic normal approach
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


ndocs = 5000
docs = gendat(nd=ndocs)

### LDA
library("topicmodels")
testLDA = posterior(LDA(docs[[1]], k=3, method = 'VEM'))  
heatmap(testLDA$topics, Rowv = NA, Colv = NA)

library(plyr); library(dplyr); library(tidyr); library(ggvis)
ddply(data.frame(treatment=docs$X$group, testLDA$topics), 'treatment', function(x) colMeans(x[,-1]) )
lm(testLDA$topics ~., data=docs$X)

gdat = gather(data.frame(docs$X, testLDA$topics), topic, value=proportion, X1, X2, X3 )

library(mgcv)
gdat %>%
  ggvis(~x, ~proportion) %>%
  group_by(topic) %>%
  layer_model_predictions(stroke=~topic, model='gam', formula = proportion ~ s(x)) 

gdat %>%
  ggvis(~group, ~proportion) %>%
  group_by(topic) %>%
  layer_model_predictions(stroke=~topic, model='lm') 
  

# testLDAdiffs = ddply(data.frame(treatment=docs$X$group, testLDA$topics), 'treatment', function(x) colMeans(x[,-1]) )  %>% 
#   select(-treatment) %>%
#   sapply(diff) %>%
#   sort %>%
#   round(3)


### STM
library(stm)
teststm0 = stm(docs[[3]], vocab=paste0('word',1:40), K=3)
teststm1 = stm(docs[[3]], vocab=paste0('word',1:40), K=3, prevalence = ~ group+x, data=docs$X)

# heatmap(teststm0$theta, Rowv = NA, Colv = NA)

ddply(data.frame(treatment=docs$X$group, teststm0$theta), 'treatment', function(x) colMeans(x[,-1]) )
ddply(data.frame(treatment=docs$X$group, teststm1$theta), 'treatment', function(x) colMeans(x[,-1]) )

lm(teststm0$theta ~., data=docs$X)
lm(teststm1$theta ~., data=docs$X)

# no model
detach(package:stm)
gdat = gather(data.frame(docs$X, teststm0$theta), topic, value=proportion, X1, X2, X3 )
gdat %>%
  ggvis(~x, ~proportion) %>%
  group_by(topic) %>%
  layer_model_predictions(stroke=~topic, model='gam', formula = proportion ~ s(x)) 

gdat %>%
  ggvis(~group, ~proportion) %>%
  group_by(topic) %>%
  layer_model_predictions(stroke=~topic, model='lm') 

# model
detach(package:stm)
gdat = gather(data.frame(docs$X, teststm1$theta), topic, value=proportion, X1, X2, X3 )
gdat %>%
  ggvis(~x, ~proportion) %>%
  group_by(topic) %>%
  layer_model_predictions(stroke=~topic, model='gam', formula = proportion ~ s(x)) 

gdat %>%
  ggvis(~group, ~proportion) %>%
  group_by(topic) %>%
  layer_model_predictions(stroke=~topic, model='lm') 


corrplot::corrplot(cor(cbind(testLDA$topics, teststm0$theta, teststm1$theta)))
ee = estimateEffect(~group, teststm1, docs$X)
plot(ee, 'group',  method = 'difference', cov.value1='control', cov.value2 = 'treatment')
ee = estimateEffect(~x, teststm1, docs$X)
plot(ee, 'x',  method = 'continuous')
