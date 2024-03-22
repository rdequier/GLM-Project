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
library(VGAM)

##### Processing Data ######

rm(list = ls())
getwd()
setwd("C:/Users/rdequ/Documents/Github/GLM-Project")

eba = read.table("eba1977.txt")


#Redefine eba as a factor
#Reference age is the group 75+
eba$age<- factor(eba$age)
eba$age <- relevel(eba$age, ref = "40-54")

#define city as a factor
eba$city <- factor(eba$city)



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



###### Exploratory Data Analysis #####

#Would be nice to update this section to make some of the graphs
#nicer using ggplot2. I am sure more can be added.

#Tables:
#There are four observations per ages group and 6 observations per city
table(eba$age)
table(eba$city)

#Frequency diagram:
plot(table(eba$cases), xlab = "cases", ylab="frequency")
#This plot shows that we can probaly rule out a zero inflated model

#Boxplots show higher cancer rate in Frederic than the other three cities
plot(eba$age, eba$cases/eba$pop)
plot(eba$city, eba$cases/eba$pop)


#Grouped histograms
ggplot(eba, aes(cases, fill=city))+
  geom_histogram(position="dodge")

#plot of observed cancer rates by city
eba_graph = eba
eba_graph$per1000 = eba$cases/eba$pop*1000

eba_graph=eba_graph %>% dplyr::select(-pop,-cases)

ggplot(eba_graph,aes(x=age, y=per1000, group=city))+
  geom_line(aes(color=city), size=0.8)+
  geom_point(aes(color=city))+
  xlab("Age Range") + ylab("Cases Per 1000") +
  ggtitle("Observed Cancer Rates by City")

#Table of average cases rates by city:
eba_graph %>% 
  dplyr::group_by(city) %>%
  dplyr::summarise(mean=mean(per1000), sd=sd(per1000))


# Note the quadratic trend in the rates, may be worth trying a quadratic model
# Note that there is a much higher rate of cases in Fredericia comapred 
# to the other two

# In the original article where this data set is used, the author mentions
# that the goal of the study was to compare one city to the other three, 
# suggesting that it might be a good idea to combine the three cities
# into one group and then 


#### Model with interaction #####
fit0 = glm(cases~city+age+city*age+offset(log(pop)),
           family=poisson(link="log"),
           data=eba)

summary(fit0)


#LR test
Anova(fit0, test = "LR",type=3)

#ineraction and age both not significant


#### Basic Model With no Interaction####              
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


#Notes: dispersion does not seem to be a big issue.
#City is not signifcant according to the wald and LRT tests


#### Basic model with continuous age #####

#use lower boundary of age
eba$NumAge=c(rep(40,4), rep(55,4), rep(60,4), rep(65,4),rep(70,4), rep(75,4))

fit2 = glm(cases~city+NumAge+offset(log(pop)),
           family=poisson(link="log"),
           data=eba)

summary(fit2) #higher devaince than before
glm.RR(fit2)

#Note the quadratic shape of the residuals in this model: justification
#For a quadratic model?
resid = resid(fit2, type="deviance")
ggplot(eba, aes(NumAge, resid)) +
  geom_point() +
  ylab("Deviance residuals") + xlab("Age")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)


resid = resid(fit2, type="pearson")
ggplot(eba, aes(NumAge, resid)) +
  geom_point() +
  ylab("Pearson residuals") + xlab("Age")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)


#### Basic Model with no city #####
fit3 = glm(cases~age+offset(log(pop)),
           family=poisson(link="log"),
           data=eba)

#### Continuous Age no City #####
fit4 = glm(cases~NumAge+offset(log(pop)),
           family=poisson(link="log"),
           data=eba)


#### Quadratic Model with continuous age #####
fit5 = glm(cases~city+NumAge+I(NumAge^2)+offset(log(pop)),
           family=poisson(link="log"),
           data=eba)

summary(fit5)



#### Quadratic Model: Continuous with no city ####
fit6=glm(cases~NumAge+I(NumAge^2)+offset(log(pop)),
         family=poisson(link="log"),
         data=eba)


#### Fredericia vs other cities: Justification ######

# Is this plot needed? Maybe we can justify using it by saying:
# We have a model that fits well according to GOF tests, but the predicted
#response shows Fredericia is much higher than the other models. This combined
#with the goal of the original dataset, which is to compare Fred to the other
#three might be a good justification.
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




#### Fredericia model #####

#Create a new column in the datafame:
eba$Fredereica = ifelse(eba$city=="Fredericia", "Fredericia", "Other")

fit7 <- glm(cases~Fredereica+age+offset(log(pop)),
            family=poisson, data=eba)
summary(fit7)




#### Fredericia Model: cts age with qaudratic #####
fit8 <- glm(cases~Fredereica+NumAge+I(NumAge^2)+offset(log(pop)),
            family=poisson, data=eba)




#### Model Selection #####

models = list(fit0,fit1, fit2, fit3, fit4, fit5,fit6,fit7, fit8)
mod.names = c("interaction", "age+city", "cts_age_city", "age_only", "age_only_cts",
              "quadratic_cts","age_only_quadratic","fredericia", "fredericia_quadratic")

aictab(cand.set=models, modnames=mod.names, second.ord=FALSE)





####### Unsure what order to do the below in:
#### Residual Plots: Basic Model #####

### Needs to be updated to apply to the models we select

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

### Needs to modified to check our models 


sim.fit1 = simulateResiduals(fit1, plot=T)

#histogram should be uniform
hist(sim.fit1)


## GOF for uniformity of the residuals and overdispersion
#This one tests uniformity of dharma residuals
testUniformity(sim.fit1)

#Compares the variance of the pearson residuals
# of the simulated vs the actual pearson residuals
testDispersion(sim.fit1)








#### Rootograms #######

#Need to update for chosen fits

#### Quasi-Likelihood to check over dispersion #####

### Needed???

#### Negative binomial model #####

#Is this needed? Need to check for overdispersion in the chosen models
fit1.nb=glm.nb(cases ~ city + NumAge +I(NumAge^2) + offset(log(pop)), data = eba)
summary(fit1.nb)



#this model does not converge:
#My guess is because the estimate of of dispersion is 1/theta, therefore
#becuase there is no quadratic trend, theta is getting a huge values





#### plotting the mean vs the variance #####

#Not sure if this is right or even necessary. Depends on whether we 
#find overdispersion in our chosen models...
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






