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
library(tidyr)


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

#deviance test function
deviance.test = function(model){
  dev=deviance(model)
  res.df = fit1$df.residual
  return(data.frame(Dev=dev, df=res.df, pvalue=(1-pchisq(dev, res.df))))
}


##### Processing Data ######


getwd()
setwd("C:/Users/rdequ/Documents/Github/GLM-Project")

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


#Grouped histograms
ggplot(eba, aes(cases, fill=city))+
  geom_histogram(position="dodge")


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

#Overdispersion
dispersiontest(fit1)

#estimating dispersion
sum(residuals(fit1, type="pearson")^2)/(fit1$df.residual)



#### Residual Plots #####

#Pearson residuals vs fitted values
resid = resid(fit1, type="pearson")
fitted=fitted(fit1)
ggplot(eba, aes(fitted, resid)) +
  geom_point() +
  ylab("Pearsons residuals") + xlab("Fitted values")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)

#Pearson residuals vs linear predictors
resid = resid(fit1, type="pearson")
linpred=fit1$linear.predictors
ggplot(eba, aes(linpred, resid)) +
  geom_point() +
  ylab("Pearsons residuals") + xlab("Linear Predictors")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)

#Deviance residuals vs fitted values#
resid = resid(fit1, type="deviance")
fitted=fitted(fit1)
ggplot(eba, aes(fitted, resid)) +
  geom_point() +
  ylab("Deviance residuals") + xlab("Fitted values")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)

#Pearson residuals vs linear predictors
resid = resid(fit1, type="deviance")
linpred=fit1$linear.predictors
ggplot(eba, aes(linpred, resid)) +
  geom_point() +
  ylab("Deviance residuals") + xlab("Linear Predictors")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)


#influentail observations:
influencePlot(fit1)


#Unsure what else we can add to this



#### Plotting the 
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

#### Fredericia vs other cities: Justification ######
summary(fit1)
beta = coef(fit1)
#Plot lines of cancer rates by city:
#Fredericia:
F = c(exp(beta[1]),
      exp(beta[1]+beta[5]),
      exp(beta[1]+beta[6]),
      exp(beta[1]+beta[7]),
      exp(beta[1]+beta[8]),
      exp(beta[1]+beta[9]))

#Horsens
H = c(exp(beta[1]+beta[2]),
      exp(beta[1]+beta[2]+beta[5]),
      exp(beta[1]+beta[2]+beta[6]),
      exp(beta[1]+beta[2]+beta[7]),
      exp(beta[1]+beta[2]+beta[8]),
      exp(beta[1]+beta[2]+beta[9]))

#Kolding
K = c(exp(beta[1]+beta[3]),
      exp(beta[1]+beta[3]+beta[5]),
      exp(beta[1]+beta[3]+beta[6]),
      exp(beta[1]+beta[3]+beta[7]),
      exp(beta[1]+beta[3]+beta[8]),
      exp(beta[1]+beta[3]+beta[9]))

#Vejle

V = c(exp(beta[1]+beta[4]),
      exp(beta[1]+beta[4]+beta[5]),
      exp(beta[1]+beta[4]+beta[6]),
      exp(beta[1]+beta[4]+beta[7]),
      exp(beta[1]+beta[4]+beta[8]),
      exp(beta[1]+beta[4]+beta[9]))

minmax = range(c(F,H,K,V))

graph_data = data.frame(x=1:6, Fredericia=F, Horsens=H, Kolden=K, Vejle = V)

graph_data = graph_data %>% pivot_longer(cols=c("Fredericia", "Horsens",
                                                "Kolden", "Vejle"),
                                         names_to="city",
                                         values_to = "rates")

labels = c("[45-54]","[55,59]","[60,64]","[65,69]","[70,74]","75+")

ggplot(graph_data, aes(x=x, y=rates))+
  geom_line(aes(color=city), size=1)+
  geom_point(aes(color=city))+
  scale_color_manual(name='City', labels=c("Fredericia", "Horsens",
                                            "Kolden", "Vejle"),
                     values=c('red', 'purple', 'steelblue','green'))+
  scale_x_continuous(breaks=1:6, labels=labels)+
  xlab("Age Range") + ylab("Predicted Cancer Rates") +
  ggtitle("Cancer Rates by City")


#Our model shows that Frederic is higher than the other three cities,
# but only Kolden is significant in the model

#Create a new column in the datafame:
eba$Fredereica = ifelse(eba$city=="Fredericia", "Fredericia", "Other")

#### Fredericia model #####
fit2 <- glm(cases~Fredereica+age+offset(log(pop)),
            family=poisson, data=eba)
summary(fit2)

#comparing the two models
anova(fit1,fit2, test="Chisq")

####Gof Tests
#pearson:
pearson.test(fit2)

#Deviance Test:
deviance.test(fit2)

#Holy trinty
#LR test
Anova(fit2, test = "LR",type=3)

#Wald
Anova(fit2, test = "Wald",type=3)

#dispersion estimate
sum(residuals(fit2, type="pearson")^2)/(fit1$df.residual)


#### plotting the mean vs the variance #####

#Not sure if this is right or even necessary#
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







#### Rootograms #######
#rootogram
rootogram(fit1,
          ylab="Root Square of Frequency",
          main = "Poisson")

rootogram(fit1.nb)

rootogram(fit2)
#### Model Selection #####

models = list(fit1, fit2)
mod.names = c("basic", "fredericia_vs_other")
aictab(cand.set=models, modnames=mod.names)




