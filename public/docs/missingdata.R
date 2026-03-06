# 14 March 2016
rm(list = ls())

library(Amelia)
library(texreg)
library(MASS)

# example data democracy on trade policy
data("freetrade") 

# Check the data to see that there is missingness on several variables
summary(freetrade)

# estimate the model from Milner and Kubota (2005)
m1 <- lm(tariff ~ polity + pop + gdp.pc + year + country, data = freetrade)
screenreg(m1)

###################################################################################
# 1st imputation
###################################################################################
# Note: include collinear variables such as different measurements of wealth (like gdp, gdp/per captia etc.)
# m = number of imputation data sets (the larger the missingness the larger the m)
# ts = time series identifier
# cs = cross section identifier
# p2s = print to screen (0 = nothing, 1 = partial, 2 = extensive) <- leave this at 2 as you can spot convergence instability easier
  # output: gives you the chain indicies before parenthesis (the last one is the chain length)
  # output: in parenthesis, the number of parameters that changed singificantly since last iteration
  # output: ! indicates a non-invertible covariance matrix
  # output: * indicates that likelihood has not monotonically increased in that step
    # many ! or * are an indication of instability
a1 <- amelia(freetrade, m = 10, ts = "year", cs = "country", p2s = 2)

# check where the imputations are stored
names(a1)

# access first imputation data set  as
# 1)
a1$imputations[[1]]
# or 2)
a1[[1]][[1]]

# quick check where the missings are
try(dev.off())
par(mfrow = c(1,1))
missmap(a1)

###################################################################################
# Combining Results: 1st) simulation; 2nd) algebra 
###################################################################################

## Simulation
# loop over imputation data sets
for (imp in 1: length(a1[[1]])){
  
  # estimate model from current imputation
  model <- lm(tariff ~ polity + pop + gdp.pc + year + country, data = a1[[1]][[imp]])
  
  # generate sampling distribution of coefficients
  S <- mvrnorm(1000, coef(model), vcov(model))
  
  # combining the simulations to incorporate imputation uncertainty
  ifelse (imp == 1, yes = S.combined <- S,
          no = S.combined <- rbind(S.combined, S))
} # end of loop over imputation data sets

# coefficients and confidence intervals 50%, 90% CI, 95% CI, 99% CI
outs <- apply(S.combined, 2, quantile, probs = c(.5, .05, .95, .025, .975, .01, 0.99))
########################################################
# result from simulation !!!
outs
########################################################

# 2nd: Task for you :)

###################################################################################
# Transformations, declaring IDs, nominal and ordinal
###################################################################################

keep.ord.scale <- "polity"
keep.nom.scale <- "signed" # 1 if country singed an IMF agreement that year
a2 <- amelia(freetrade, m = 10, ts = "year", cs = "country", p2s = 2,
             ords = keep.ord.scale, noms= keep.nom.scale)

# see that the ordinal scale of the polity variable is kept in a2 (this is not recommended)
table(a1[[1]][[1]]$polity)
table(a2[[1]][[1]]$polity)

# see that the nominal scale of the signed variable is kept: nominal variables must be declared (e.g. denomination, ethnicity)
# also do this if you need your DV to be binary
table(a1[[1]][[1]]$signed)
table(a2[[1]][[1]]$signed)

## heavily skewed variables should be log tansformed
try(dev.off())
par(mfrow = c(1,2))
hist(freetrade$tariff)
hist(log(freetrade$tariff))

hist(freetrade$pop)
hist(log(freetrade$pop))

hist(freetrade$gdp.pc)
hist(log(freetrade$gdp.pc))

## transformations that can be set in Amelia
# sqrts = take the square root for count data
# log = take the natural log for heavily skewed data
# lgstc = take a logistic transformtion for proportions

log.trans <- c('tariff', 'pop', 'gdp.pc')
a3 <- amelia(freetrade, m = 10, ts = "year", cs = "country", p2s = 2,
             ords = keep.ord.scale, noms= keep.nom.scale,
             logs = log.trans)

# declare ID variables or any variables that you do not want to include in the imputation model
# if you do not declare character variables, Amelia will return an error
freetrade$new.id <- seq(from = 1, to = nrow(freetrade), by = 1)
freetrade$coders.comments <- rep("useful comments by the coder", nrow(freetrade))

declare.ids <- c('new.id', 'coders.comments')
a4 <- amelia(freetrade, m = 10, ts = "year", cs = "country", p2s = 2,
             ords = keep.ord.scale, noms= keep.nom.scale, logs = log.trans,
             idvars = declare.ids)

###################################################################################
# Time Series & Panel
###################################################################################
# if tariffs vary smoothly over time, some polynomial exists that describes tariffs in unit i at time t:
# tariff_i_j = beta0 + beta1 *t + beta2 *t^2 + beta3 *t^3 ....
a5 <- amelia(freetrade, m = 10, ts = 'year', cs = 'country', p2s = 2, 
             idvars = declare.ids, logs = log.trans,
             polytime = 2)

# to allow the time effects to vary over country we intersect cross-sectional unit with the polynomials
a6 <- amelia(freetrade, m = 10, ts = 'year', cs = 'country', p2s = 2, 
             idvars = declare.ids, logs = log.trans, polytime = 2,
             intercs = TRUE)

# we can run fixed effects by setting intersecs without polytime
a7 <- amelia(freetrade, m = 10, ts = 'year', cs = 'country', p2s = 2, 
             idvars = declare.ids, logs = log.trans,
             intercs = TRUE)

# lags and leads
a8 <- amelia(freetrade, m = 10, ts = 'year', cs = 'country', p2s = 2, 
             idvars = declare.ids, logs = log.trans, intercs = TRUE,
             lags = 'tariff', leads = 'tariff')

## evaluate the predictions in our time series
dev.off()
par(mfrow = c(2,2))
tscsPlot(a4, cs = 'Malaysia', main = 'Malysia, no time settings', var = 'tariff', ylim = c(0, 50))
tscsPlot(a5, cs = 'Malaysia', main = 'Malysia, time effect', var = 'tariff', ylim = c(0, 50))
tscsPlot(a6, cs = 'Malaysia', main = 'Malysia, country specific time effect', var = 'tariff', ylim = c(0, 50))
tscsPlot(a8, cs = 'Malaysia', main = 'Malysia, Lags and Leads', var = 'tariff', ylim = c(0, 50))

# use splines
a9 <- amelia(freetrade, m = 10, ts = 'year', cs = 'country', p2s = 2, 
             idvars = declare.ids, logs = log.trans, intercs = TRUE,
             splinetime = 2)
# polynomials vs splines (splines are much faster, require less parameters to be estimated)
dev.off()
par(mfrow = c(1,2))
tscsPlot(a6, cs = 'Malaysia', main = 'Malysia, time polynomials X country', var = 'tariff', ylim = c(0,30))
tscsPlot(a9, cs = 'Malaysia', main = 'Malysia, splines X country', var = 'tariff', ylim = c(0,30))

###################################################################################
# Priors
###################################################################################
# use ridge priors for high degree of missingness or strong correlations among variables or
# number of observations is only slightly larger than number of variables

# a ridge prior adds numerical stability by shrinking covariances among variables towards 0
# this is roughtly equivalent with adding x new observations with same means and variances but 0 covariance
# reduces variance but increases bias
# a ridge prior of 1% is a starting point, 5% is okay, 10% is the upper bound

# estimate time series and cross-section model with 5% ridge prior
# see how chainlengths decrease rapidly
a10 <- amelia(freetrade, m = 10, ts = 'year', cs = 'country', p2s = 2, idvars = declare.ids, 
              polytime = 2, intercs = TRUE,
              empri = .01 * nrow(freetrade))
par(mfrow= c(1,2))
tscsPlot(a6, cs =  'Malaysia', main = 'Malaysia without ridge prior', var = 'tariff', ylim = c(0,30))
tscsPlot(a10, cs =  'Malaysia', main = 'Malaysia with ridge prior', var = 'tariff', ylim = c(0,30))

## Observation-level priors
freetrade[freetrade$country=='Thailand', c('year', 'country', 'tariff')]

# we have expert knowlege: 
# tariff rates were roughly 40% in 1986 to 1988
# the expert says: plus minus 3%

## Option 1
# columns: 1st = row numbers, 2nd = column number of variable in data, 3rd = mean, 4th = standard deviation
prior.mat <- matrix(c(158, 159, 160, 3, 3, 3, 40, 40, 40, 3, 3, 3), nrow = 3, ncol = 4)

## Option 2
# columns: 1st = row numbers, 2nd = column number of variable in data, 3rd = lower bound, 4th = upper bound,
# 5th = confidence level
prior.mat2 <- matrix(c(158, 159, 160, 3, 3, 3, 34, 34, 34, 46, 46, 46, .95, .95, .95),
                     nrow = 3, ncol = 5)

# estimate imputation with expert priors
# if you keep the log transformation option, you need to log transform the priors as well!!
a11 <- amelia(freetrade, m = 10, ts = 'year', cs = 'country', p2s = 2, 
              idvars = declare.ids, polytime = 2, intercs = TRUE,
              priors = prior.mat)
tscsPlot(a11, cs =  'Thailand', main = 'Thailand with expert priors', var = 'tariff',
         ylim = c(-10, 60))

###################################################################################
# Diagnostics
###################################################################################

## Overlay Densities
# Note that densities may look different as our MAR assumption states that missingness
# is random conditional on the data in the imputation model
# However, when you see big differences think whether it is plausible that certain values 
# are more likely missing than others
dev.off()
par(mfrow = c(1,3))
compare.density(a6, var = "tariff")
compare.density(a6, var = "polity")
compare.density(a6, var = "signed")

# in time series you can use the tscsPlot we used before but the same note of caution because of MAR applies

## Overimpute
# In each step 1 observed value in the data is dropped and then predicted using all other data to predict
# (this is done several hundred times to more accuratley predict confidence intervals)
# If you are interested overimputation can be used to get at measurement error
dev.off()
par(mfrow = c(1,1))
# the intervals plotted are 90% CI's
# the colors correspond to the percentage missing on that observation
overimpute(a6, var = 'polity')

## Overdispersed Starting Values
# EM may converge on a local maximum - to check you can choose to use a number of 
# different starting values and see if the alogrithm converges on the same maximum
# that was found in your intial run.
# A lot of parameters are imputed (high dimensionality), therefore you can choose to
# display the larges or the do largest principal components (dims)

# larges principal component; 6 new starting values
disperse(a6, dims = 1, m = 6)

# 2 largest pc's; 6 new starting values
disperse(a6, dims = 2, m = 6)

# multiple points of conversion indicate your model is not identified 
freetrade2 <- freetrade
freetrade2$tariff2 <- freetrade2$tariff * 2 + 3 
a12 <- amelia(freetrade2, m = 10, ts = 'year', cs = 'country', p2s = 2,
              idvars = declare.ids)
disperse(a12, dims = 2, m = 6)
