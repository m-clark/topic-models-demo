library(logitnorm)

x = rlogitnorm(n=10000)  # mu 0 sd 1
head(x)
ggplot2::qplot(x, geom='density')
ggplot2::qplot(qlogis(x), geom='density')

x = rlogitnorm(n=10000, sigma=3)  # mu 0 sd 1
head(x)
ggplot2::qplot(x, geom='density')
ggplot2::qplot(qlogis(x), geom='density')


# note that rlogitnorm is just plogis(rnorm)

# multivariate
library(MASS)
cormat = matrix(c(1,.5,.5,
                  .5,1,.1,
                  .5,.1,1), 3, 3)
x = mvrnorm(1000, mu=rep(0,3), Sigma=cormat, empirical = T)

x_prob = t(apply(x, 1, function(row) exp(row)/sum(exp(row))))
head(x_prob)
head(gtools::rdirichlet(100, rep(.5,3)))


X = cbind(1, rnorm(100), rnorm(100))

w1 = c(0, 1, -1)
w2 = c(0, -1, 1)
w3 = c(0, .1, .1)

topicProb0 = cbind(X%*%w1, X%*%w2, X%*%w3)

logisticNormal = t(apply(topicProb0, 1, function(row) plogis(row)/sum(plogis(row))))
soft = t(apply(topicProb0, 1, function(row) exp(row)/sum(exp(row))))

comparison = data.frame(ln=logisticNormal, softmax=soft)
head(comparison)
corrplot::corrplot(cor(comparison), method='number')
