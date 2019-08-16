rm(list = ls())
setwd("C:\\Users\\phili\\Dropbox\\Intro to QM 2015 SMLL\\Data")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(foreign)
df <- read.dta("MROZ.dta")

# generate variable kids (either have kids or not)
df$kids <- NA
df$kids[df$kidslt6 > 0 | df$kidsge6 > 0] <- 1
df$kids[df$kidslt6 == 0 & df$kidsge6 == 0] <- 0
table(df$kids)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Logistic regression
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model
m1 <- glm(inlf ~ kids + age + educ, data=df, family=binomial(logit))
summary(m1)

# make a prediction for: woman, no kids, 32 yrs. old, 13 yrs. of education
coef(m1) # all coeficients
coef(m1)[1] # first (i.e. intercept)
coef(m1)[2]
# latent y
y_hat <- coef(m1)[1] + 0 * coef(m1)[2] + coef(m1)[3] * 32 + coef(m1)[4] * 13
names(y_hat) <- "y_hat"
y_hat
# put y_hat into the link function
1 / ( 1 + exp(- y_hat) ) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Interpretation 1: Prediction for 1 ideal type (avg. woman)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Zelig)

# step 1: estimate model
z.out1 <- zelig(inlf ~ kids + age + educ + exper + huseduc + huswage, data = df,
                model = "logit", cite = FALSE)

# step 2: set covariate values (avg. woman)
avg.woman <- setx(z.out1, kids = median(df$kids), age = mean(df$age), educ = mean(df$educ),
                  exper = mean(df$exper), huseduc = mean(df$huseduc), huswage = mean(df$huswage))

# step 3: simulate
sim.out1 <- sim(z.out1, avg.woman)
summary(sim.out1)
plot(sim.out1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Interpretation 2: Compare 2 groups (kids and no kids)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model already estimated, no need to do that again

# set covariates
no.kids <- setx(z.out1, kids = 0, age = mean(df$age), educ = mean(df$educ),
                exper = mean(df$exper), huseduc = mean(df$huseduc), huswage = mean(df$huswage))
has.kids <- setx(z.out1, kids = 1, age = mean(df$age), educ = mean(df$educ),
                 exper = mean(df$exper), huseduc = mean(df$huseduc), huswage = mean(df$huswage))

# simulate
sim.out2 <- sim(z.out1, no.kids, has.kids)
summary(sim.out2) # Check the first differences!!!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Interpretation 3: Show effect of a continuous variable (education - range: 5 to 17)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(df$educ)
# set covariates
set.seed(123)
edu.range <- setx(z.out1, kids = median(df$kids), age = mean(df$age), educ = 5:17,
                  exper = mean(df$exper), huseduc = mean(df$huseduc), huswage = mean(df$huswage))
# simulate
sim.out3 <- sim(z.out1, edu.range)
par(mfrow = c(1,1))
ci.plot(sim.out3, qi =  "ev", ci = 95, main = "Effect of Education on Employment")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model Quality
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
predicted.probabilities <- predict(m1, type = "response")
expected.values <- ifelse(predicted.probabilities > 0.5, yes = 1, no = 0)
observed.outcomes <- m1$model$inlf # that's the dependent variable

# prediction quality table
out.table <- table(observed.outcomes, expected.values)
out.table

# correctly predicted
corr.pred <- ( sum(diag(out.table)) / sum(out.table) )
corr.pred

# the naive guess
summary(df$inlf)
# the modal category is 1
# the mean is .57
# therefore, always saying a women works will make you right 57% of the time


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Joint Hypothesis Testing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# first model
m1 <- glm(inlf ~ kids + age + educ + exper, data=df, family=binomial(logit))
m2 <- glm(inlf ~ kids + age + educ + exper + huseduc + huswage, data=df, family=binomial(logit))
library(texreg)
screenreg(list( m1, m2))

# notice that houseduc and huswage are insignificant, but are they jointly significant
anova(m1, m2, test = "Chisq")

# alternatively, you can go by AIC and BIC
# lower values are better
# BIC is more conservative
AIC(m1, m2)
BIC(m1, m2)
