X = cbind(1, rnorm(100), rnorm(100))

w1 = c(0, 1, -1)
w2 = c(0, -1, 1)
w3 = c(0, .1, .1)

topicProb0 = cbind(X%*%w1, X%*%w2, X%*%w3)

logisticNormal = t(apply(topicProb0, 1, function(row) plogis(row)/sum(plogis(row))))
soft = t(apply(topicProb0, 1, function(row) exp(row)/sum(exp(row))))

head(cbind(logisticNormal, soft))
corrplot::corrplot(cor(cbind(logisticNormal, soft)), method='number')
