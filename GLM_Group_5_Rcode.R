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
library(ggpubr)
library(kableExtra)
library(ggeffects)
library(see)
library(topmodels)
library(vcd)
library(performance)
##### Processing Data ######

rm(list = ls())
getwd()
setwd("C:/Users/rdequ/Documents/Github/GLM-Project")

eba = read.table("eba1977.txt")


#Redefine eba as a factor
#Reference age is youngest group
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
  res.df =model$df.residual
  return(data.frame(Dev=dev, df=res.df, pvalue=(1-pchisq(dev, res.df))))
}



###### Exploratory Data Analysis #####

#Would be nice to update this section to make some of the graphs
#nicer using ggplot2. I am sure more can be added.


#Does the model follow a poisson, binomial, or negative binomial?
distplot(eba$cases, type="poisson")
distplot(eba$cases, "nbinomial")
distplot(eba$cases, "binomial")
#There is evidence for both poisson and negative binomial


#Check mean vs variances of data
mean(eba$cases)
var(eba$cases) #very close: overdispersion may not be an issue

sqrt(mean(eba$cases))
sd(eba$cases)


#Cases by age and group
eba %>% group_by(city) %>% 
  dplyr::summarise(mean=mean(cases),
                   var=var(cases),
                   sd=sd(cases),
                   med=median(cases))

eba %>% group_by(age) %>% 
  dplyr::summarise(mean=mean(cases),
                   var=var(cases),
                   sd=sd(cases),
                   med=median(cases))



# We are more interested in comparing rates
#Add new columns:
eba$rate = eba$cases/eba$pop
eba$perthousand = eba$rate*1000

rates_city = eba %>% group_by(city) %>% 
  dplyr::summarise(mean=mean(perthousand),
                   median=median(perthousand),
                   sd=sd(perthousand),
                   var=var(perthousand)) %>% 
  mutate_if(is.numeric, round,digits=3) 
rates_city

#export to latex
kable(rates_city,"latex", booktabs=TRUE)


#Comparing Rates vs Age
rates_age = eba %>% group_by(age) %>% 
  dplyr::summarise(mean=mean(perthousand),
                   median=median(perthousand),
                   sd = sd(perthousand),
                   var=var(perthousand)) %>% 
  mutate_if(is.numeric, round,digits=3) 
rates_age
kable(rates_age, "latex", booktabs=TRUE)


#Tables:
#There are four observations per ages group and 6 observations per city
table(eba$age)
table(eba$city) 

#Frequency diagram:
plot(table(eba$cases), xlab = "cases", ylab="frequency")
#This plot shows that we can probably rule out a zero inflated model

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

#goodness of fit tests show this model does not fit well
pearson.test(fit0)
deviance.test(fit0)


#LR test
Anova(fit0, test = "LR",type=3)

#There is no evidence of an interaction


#### Basic Model With no Interaction####              
#fit a basic model with offset for population sizes
fit1 = glm(cases~city+age+offset(log(pop)),
               family=poisson(link="log"),
               data=eba)

summary(fit1)


#Gof Tests
#pearson:
pearson.test(fit1)

#Deviance Test:
deviance.test(fit1)

#LR test
Anova(fit1, test = "LR",type=3)

#Wald
Anova(fit1, test = "Wald",type=3)


#estimating dispersion
sum(residuals(fit1, type="pearson")^2)/(fit1$df.residual)


#Notes: dispersion does not seem to be a big issue.
#City is not significant according to the wald and LRT tests
#However, in the model summary, we see that
#We also note a dispersion parameter of 1.5


# Other possibilities: negative binomial
fit1.nb=glm.nb(cases~city+age+offset(log(pop)),
               data=eba)
#Note that theta diverges: indicating there is no quadratic
#trend in the variance


#There are no zeros in the model: try a zero-truncated distribution
fit1.zero=zerotrunc(cases~city+age+offset(log(pop)),
                    dist="poisson",
                    data=eba)
summary(fit1.zero)
AIC(fit1.zero)
#we have a nearly identical estimates and AIC. Does not improve fit, but
#complicated model interpretation.

#### Basic model with continuous age #####
#Age groups are unevenly spaced, so might be worth trying:
#use lower boundary of age because the highest age group
eba$NumAge=c(rep(40,4), rep(55,4), rep(60,4), rep(65,4),rep(70,4), rep(75,4))

fit2 = glm(cases~city+NumAge+offset(log(pop)),
           family=poisson(link="log"),
           data=eba)

summary(fit2) #higher deviance than before

pearson.test(fit2)
deviance.test(fit2)
#Worse fit, but Note the quadratic shape of the residuals in this model: 
# justification for a quadratic model?
resid = resid(fit2, type="deviance")
r1=ggplot(eba, aes(NumAge, resid)) +
  geom_point() +
  ylab("Deviance residuals") + xlab("Age")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)+
  ggtitle("Deviance Residuals")+
  theme(plot.title = element_text(hjust = 0.5))


resid = resid(fit2, type="pearson")
r2=ggplot(eba, aes(NumAge, resid)) +
  geom_point() +
  ylab("Pearson residuals") + xlab("Age")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)+
  ggtitle("Pearson Residuals")+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(r1,r2)




#### Factor Model with no city #####
#Since city was not significant:
fit3 = glm(cases~age+offset(log(pop)),
           family=poisson(link="log"),
           data=eba)
summary(fit3)
pearson.test(fit3)
deviance.test(fit3)
#This model is borderline with deviance goodness of fit test

#### Continuous Age no City #####
fit4 = glm(cases~NumAge+offset(log(pop)),
           family=poisson(link="log"),
           data=eba)

pearson.test(fit4)
deviance.test(fit4)
#Very bad fit

#### Quadratic Model with continuous age #####
fit5 = glm(cases~city+poly(NumAge,2)+offset(log(pop)),
           family=poisson(link="log"),
           data=eba)
pearson.test(fit5)
deviance.test(fit5)
#Acceptable fit

#### Quadratic Model: Continuous with no city ####
fit6=glm(cases~poly(NumAge,2)+offset(log(pop)),
         family=poisson(link="log"),
         data=eba)
pearson.test(fit5)
deviance.test(fit5)

#### Fredericia vs other cities: Justification ######

#Recall there was a significant difference between Fredericia 
#and Kolding, while the difference between Horsens and Fredericia
#was borderline signficant
summary(fit1)

#Note the very small difference between the three cities 
#compared to Fredericia
plot(predict_response(fit1, c("age","city"), condition = c(pop=100)))

#In our EDA, we also say a notably higher rate of cancer cases
#compared to the other three cities:
rates_age
rates_city

#We therefore group the three other cities and compare them to Fredericia

#### Fredericia model #####

#Create a new column in the datafame:
eba$Fredericia = ifelse(eba$city=="Fredericia", 1, 0)

fit7 <- glm(cases~Fredericia+age+offset(log(pop)),
            family=poisson, data=eba)

summary(fit7)

#Deviance is lower, and so is AIC

#### Fredericia Model: cts age with qaudratic #####
fit8 <- glm(cases~Fredericia+poly(NumAge,2)+offset(log(pop)),
            family=poisson, data=eba)

summary(fit8)

#The quadratic modle with binary city has the lowest AIC
#and 


#### Model Selection #####

#omit continuous models with age, unless age is squared:
#these are fit2 and fit4

models = list(fit1, fit3, fit5,fit6,fit7, fit8)
mod.names = c("Factor: Age+City", "Factor: Age",
              "Continuous: Age^2+City","Continuous: Age^2",
              "Factor: Binary City", "Continuous: Binary City+Age^2")

aic_table = aictab(cand.set=models, modnames=mod.names, second.ord=FALSE)
aic_table


#Two best models are the Linear Factor Model, and Quadratic Continuous Model,
#Both comparing Fredericia to the other three cities

summary(fit7)
summary(fit8)

#GoF tests for the two best models:
pearson.test(fit7)
pearson.test(fit8)

deviance.test(fit7)
deviance.test(fit8)


#### Residuals of Two Best Models: #####

#Two best models both contain only the city Fredericia:

#deviance residuals vs linear predictors for fits 7 and 8:
resid = resid(fit7, type="deviance")
linpred=fit7$linear.predictors
p1=ggplot(eba,aes(linpred, resid)) +
  geom_point() +
  ylab("Deviance residuals") + xlab("Linear Predictors")+
  ggtitle("Linear Factor: Deviance")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)

resid = resid(fit8, type="deviance")
linpred=fit8$linear.predictors
p2=ggplot(eba,aes(linpred, resid)) +
  geom_point() +
  ylab("Deviance residuals") + xlab("Linear Predictors")+
  ggtitle("Quadratic: Deviance")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)


#Cook's distance vs index for fits 7 and 8:
cooks=cooks.distance(fit7)
p3=ggplot(data.frame(x=seq(cooks),y=cooks))+
  geom_line(aes(x,y), col="red")+
  geom_point(aes(x,y)) +
  ylab("Cook's Distance") + xlab("Index")+
  ggtitle("Cook's Distance vs Index")+
  ylim(0,0.6)


cooks=cooks.distance(fit8)
p4=ggplot(data.frame(x=seq(cooks),y=cooks))+
  geom_line(aes(x,y), col="red")+
  geom_point(aes(x,y)) +
  ylab("Cook's Distance") + xlab("Index")+
  ggtitle("Cook's Distance vs. Index")


#Linear predictors vs age
resid = resid(fit7, type="deviance")
p5=ggplot(eba, aes(NumAge, resid)) +
  geom_point() +
  ylab("Deviance residuals") + xlab("Age")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)+
  ggtitle("Deviance Residuals vs Age")


resid = resid(fit8, type="deviance")
p6=ggplot(eba, aes(NumAge, resid)) +
  geom_point() +
  ylab("Deviance residuals") + xlab("Age")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)+
  ggtitle("Deviance Residuals vs Age")


ggarrange(p1,p3,p5,p2,p4,p6)

#The models perform similarly, but the quadratic model has two influential
#observations


#### DHARMa residuals to check for over dispersion ####

#Estimated dispersion parameters for both models:
sum(residuals(fit7, type="pearson")^2)/(fit7$df.residual)
sum(residuals(fit8, type="pearson")^2)/(fit8$df.residual)


#Tests for the Linear Factor Model
sim.fit7 = simulateResiduals(fit7, plot=T)
sim.fit8 = simulateResiduals(fit8, plot=T) #outlier detcted
hist(sim.fit7)
hist(sim.fit8) # outlier detected

#The models perform similarly. We choose the model, fit7, that treats
#age as a factor for two reasons: there are influential observaitons
#and outliers in the quadratic model. The factor model is also 
#much easier to interpret than a quadratic poisson model with an 
#orthognal polynomial structure



###### Diagnostics of Chosen model ######

#Observed vs expected counts and rates per 1000
topmodels::rootogram(fit7, style="standing")

fit7.predictions = cbind(eba[c("city", "age","pop","cases")],
                         Predicted=round(exp(predict(fit7)),2))
fit7.predictions



#Checking residual plots and for influential observations:
#Nothing detected
resid = resid(fit7, type="pearson")
linpred=fit7$linear.predictors
p1=ggplot(eba,aes(linpred, resid)) +
  geom_point() +
  ylab("Pearson residuals") + xlab("Linear Predictors")+
  ggtitle("Pearson Residuals vs Linear Predictors")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)


resid = resid(fit7, type="deviance")
linpred=fit8$linear.predictors
p2=ggplot(eba,aes(linpred, resid)) +
  geom_point() +
  ylab("Deviance residuals") + xlab("Linear Predictors")+
  ggtitle("Deviance Residuals vs. Liear Predictors")+
  geom_smooth()+
  geom_hline(yintercept=0, col="red", lty=2)


cooks=cooks.distance(fit7)
p3=ggplot(data.frame(x=seq(cooks),y=cooks))+
  geom_line(aes(x,y), col="red")+
  geom_point(aes(x,y)) +
  ylab("Cook's Distance") + xlab("Index")+
  ggtitle("Cook's Distance vs. Index")


ggarrange(p1,p2,p3, ncol=1)


#mean-variance relationship:
check_model(fit7, check="overdispersion")
#The mean variance relationship is modeled acceptably, but
#there are some deviations from a linear trend through the origin.


#Dispersion: 
sum(residuals(fit7, type="pearson")^2)/(fit7$df.residual)
#A value of 1.33 is a low level of dispersion, but should
#still be cocrrecteed using either the sandwich estimator or
#a quasipoisson model. We have already seen that a negative bimomial
#model does not improve the fit. 

#For more evidence of this, examine the current model's mean 
#variance relationship with a negative binomial distribution:

fit7.nb=glm.nb(formula = cases ~ Fredericia + age + offset(log(pop)),
         data = eba)
check_model(fit7.nb, check="overdispersion")




#### QuasiPoisson + Sandwich Estimator ####

#Sandwich Estimator on normal model:
coeftest(fit7, vcov=sandwich)
#Note that the sandwich estimaor underestimates standard errors.
#This is because of small sample size

#QausiPoisson:

fit7.quasi = glm(cases~Fredericia+age+offset(log(pop)),
                 family=quasipoisson, data=eba)
summary(fit7.quasi)
#coefficient test with robust estimator


#LR +  Wald Tests
Anova(fit7.quasi, test = "Wald",type=3)
Anova(fit7.quasi, test = "LR",type=3)

#We get the same estimates as with the Poisson model, but now
#the effect of city is borderline significant.

