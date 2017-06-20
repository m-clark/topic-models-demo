library(DirichletReg) # for rdirichlet

### Explicit representation, as depicted in the usual plate diagram
Ndocs = 500                                                                     # Number of documents
WordsPerDoc = rpois(Ndocs, 40)                                                  # Number of total words in a document
thetaList = list( c(.15,.25,.6), c(.1, .1, .8))                                 # Topic proportions for first and second half of data
theta1 = t(replicate(Ndocs/2, thetaList[[1]]))                                  # First half of documents show higher probability of topic '1' and '2'
theta2 = t(replicate(Ndocs/2, thetaList[[2]])) 
theta = rbind(theta1, theta2)                                                   # Topic probabilities for 500 docs

Z = t(apply(theta, 1, function(topprob) rmultinom(1, 1, topprob)))              # draw topic assignment
# colMeans(Z[1:Ndocs/2,])   # roughly equal to theta1

z = apply(Z, 1, function(row) which(row==1))                                    # topic assignment as label 1:3

B_k = rdirichlet(n=40, alpha=rep(.05,3))                                        # topics, i.e. distribution over words
# B_k 

wordlist1 = sapply(1:Ndocs, function(i) t(rmultinom(1, WordsPerDoc[i], B_k[,z[i]]))  # given topic assignment, draw words according to doc size
                  , simplify = F)  
ldadat1 = do.call(rbind, wordlist1)                                                  # smash to doc-term matrix
wordlist1 = lapply(wordlist1, function(wds) rep(paste0('word',1:length(wds)), wds))   # bag of words representation
# table(wordlist1[[1]])


### Matrix approach, don't need explicit topic label
ZB = tcrossprod(Z, B_k)                                                
ldadat2 = t(sapply(1:Ndocs, function(i) rmultinom(1, WordsPerDoc[i], ZB[i,])))
wordlist2 = apply(ldadat2, 1, function(row) rep(paste0('word', which(row!=0)), row[row!=0]) )


# curiosity, overlap of words
nWordOverlap = rep(NA, Ndocs)
percWordOverlap = rep(NA, Ndocs)
for (i in 1:Ndocs){
  nWordOverlap[i] = length(intersect(wordlist1[[i]], wordlist2[[i]]))
  percWordOverlap[i] = length(intersect(wordlist1[[i]], wordlist2[[i]])) /
                       length(unique(c(wordlist1[[i]], wordlist2[[i]])))
}
summary(nWordOverlap)
summary(percWordOverlap)


### Run models
# depending on Ndocs and other, it is possible to get an unsatisfactory/random
# result; usually a redo is enough to recover topic probs; the settings are an
# attempt to get better result; but will slow things down

library(topicmodels)
controlSettings = list(nstart=50) #, verbose=75, var=list(iter.max=5000, tol=10e-8), em=list(iter.max=5000, tol=10e-8)
# controlSettings = NULL
LDA1 =LDA(ldadat1, k=3, method = 'VEM', control=controlSettings)
LDA2 = LDA(ldadat2, k=3, method = 'VEM', control=controlSettings)
LDA1Post = posterior(LDA1)
LDA2Post = posterior(LDA2)
# heatmap(LDA1Post$topics, Rowv = NA, Colv = NA)
# heatmap(LDA2Post$topics, Rowv = NA, Colv = NA)

firsthalf = 1:(Ndocs/2)
secondhalf = (Ndocs/2+1):Ndocs

round(
  rbind(thetaList[[1]], thetaList[[2]], 
      sort(colMeans(LDA1Post$topics[firsthalf,])), sort(colMeans(LDA1Post$topics[secondhalf,])),
      sort(colMeans(LDA2Post$topics[firsthalf,])), sort(colMeans(LDA2Post$topics[secondhalf,])) )
  , 2)

library(LDAvis)
shinyJSON = createJSON(phi=exp(LDA1@beta), theta=LDA1@gamma, doc.length = WordsPerDoc, 
                       vocab = paste0('word',1:40), term.frequency = colSums(ldadat1))
serVis(shinyJSON)
