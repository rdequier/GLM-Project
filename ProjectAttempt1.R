### GLM PROJECT ###

#### Packages #####
library(dplyr)
library(ggplot2)
library(car)
library(countreg)
library(AER)
library(AICcmodavg)
library(ggplot2)
library(DHARMa)
library(topmodels)
library(MASS)
library(sandwich)
library(lmtest)

######Functions########

#function to create relative risk confidence intervals
glm.RR = function(model){
  
  print(exp(cbind(RR = model$coefficients, confint(model))),
        digts=3)
}

#Function to cmopute pearson goodnoess of fit test
pearson.test= function(model){
  X2 = sum(residuals(model, type="pearson")^2)
  n = dim(eba)[1]
  p = length(coef(model))
  return(data.frame(X2s=X2, pvalue=(1-pchisq(X2,n-p))))
}

deviance.test = function(model){
  dev=deviance(model)
  res.df = fit1$df.residual
  return(data.frame(Dev=dev, df=res.df, pvalue=(1-pchisq(dev, res.df))))
}


#Function to plot deviance residuals vs linear predictors:
linpred_vs_dev = function(model){
model.rdev<-residuals(model,type="deviance")
plot(predict(model),model.rdev,xlab="Linear predictor", ylab="Devianceresiduals")
abline(h=0,lty=2,col="grey")
loess.model.rdev <- loess(model.rdev~predict(model))
model.lo.pred.dev <- predict(loess.model.rdev, se=T)
j.model<- order(predict(model))
lines(predict(model)[j.model],model.lo.pred.dev$fit[j.model],col="blue",lwd=3)
lines(predict(model)[j.model],model.lo.pred.dev$fit[j.model]
      +2*model.lo.pred.dev$s[j.model], lty=2,col="red")
lines(predict(model)[j.model],model.lo.pred.dev$fit[j.model]
      -2*model.lo.pred.dev$s[j.model], lty=2,col="red")
}



##### Processing Data ######


getwd()
setwd("C:/Users/rdequ/Documents/R/GLM")

eba = read.table("eba1977.txt")


#Redefine eba as a factor
#Reference age is the group 75+
eba$age<- factor(eba$age)
eba$age <- relevel(eba$age, ref = "40-54")

#define city as a factor
eba$city <- factor(eba$city)


###### Exploratory Data Analysis #####


#Tables:
#There are four observations per ages group and 6 observations per city
table(eba$age)
table(eba$city)

#Frequency diagram:
plot(table(eba$cases), xlab = "cases", ylab="frequency")

#Boxplots show higher cancer rate in Frederic than the other three cities
plot(eba$age, eba$cases/eba$pop)
plot(eba$city, eba$cases/eba$pop)





#### Basic Model With Offset ####              
#fit a basic model with offset for population sizes
fit1 = glm(cases~city+age+offset(log(pop)),
               family=poisson(link="log"),
               data=eba)

summary(fit1)
glm.RR(fit1)


####Gof Tests
#pearson:
pearson.test(fit1)

#Deviance Test:
deviance.test(fit1)

#Holy trinty
#LR test
Anova(fit1, test = "LR",type=3)

#Wald
Anova(fit1, test = "Wald",type=3)

#Deviance residuals vs linear predictor:
linpred_vs_dev(fit1)


#Overdispersion
dispersiontest(fit1, alternative ="less")
dispersiontest(fit1)




#### DHARMa residuals to check for over dispersion ####
sim.fit1 = simulateResiduals(fit1, plot=T)

#histogram should be uniform
hist(sim.fit1)


## GOF for uniformity of the residuals and overdispersion
#This one tests uniformity of dharma residuals
testUniformity(sim.fit1)

#Compares the variance of the pearson residuals
# of the simulated vs the actual pearson residuals
testDispersion(sim.fit1)


#### Quasi-Likelihood to check over dispersion #####
fit1.quasi =  glm(cases ~ city + age, offset = log(pop),
                  family = quasipoisson(link = "log"), data = eba)
summary(fit1.quasi)

#Dispersion parameter is estimated to be 1.5
#Notice that if we compare the summaries for the two models,
#the significance changes
summary(fit1.quasi)
summary(fit1)

#With the quasi-likelihood model, city-Kolding is no longer significant,
#so we may be at the risk of type I error

#### Negative binomial model #####
fit1.nb=glm.nb(cases ~ city + age + offset(log(pop)), data = eba)
summary(fit1.nb)
#this model does not converge:
#My guess is because the estimate of of dispersion is 1/theta, therefore
#becuase there is no quadratic trend, theta is getting a huge values



#### Basic fit with robust estimator ####
coeftest(fit1, vcov=sandwich)

#Now the cities of Kolding and Vejle are significant.. hmmm

#### plotting the mean vs the variance #####

#Not sure if this is right....?
xb <- predict(fit1.nb)
g <- cut(xb, breaks=quantile(xb,seq(0,100,10)/100))
m <- tapply(eba$cases, g, mean)
v <- tapply(eba$cases, g, var)
plot(m, v, xlab="Mean", ylab="Variance", main="Mean-Variance Relationship")
mtext("Articles Published by Ph.D. Biochemists",padj=-0.5)
x <- seq(6,10,0.02)
lines(x, x*summary(fit1.quasi)$dispersion, lty="dashed")
lines(x, x*(1+x/fit1.nb$theta))
legend("topleft", lty=c("dashed","solid"),legend=c("Q. Poisson","Neg. Binom."), inset=0.05)




#### Model with Fredericia as the only city ####
fit2 = glm(cases~1+offset(log(pop)),
           family=poisson(link="log"),
           data=eba)

summary(fit2)


#Not sure how to do this.


#### Rootograms #######
#rootogram
rootogram(fit1,
          ylab="Root Square of Frequency",
          main = "Poisson")

rootogram(fit1.nb)
#### Model Selection #####

models = list(fit1, fit1.nb)
mod.names = c("basic", "nb")
aictab(cand.set=models, modnames=mod.names)




