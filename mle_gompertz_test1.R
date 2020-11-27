## mle_gompertz_test.R

## we have seen the MLE function fit simulated Gompertzian data well, so the next step is to try it out on observed data
## Notes:
## * in the simulation, we assigned "true" values for the parameters M and beta, which are theoretically unobservable
## * for this test, parameter values will be replaced by "guesses" for M and beta based on what the observed numident data tell us (likewise, u and l will be observed min and max numident death_age values)
## * whereas in the simulation ages at death were on a continuous real-number range, here numident death_age takes on floored integer values. This discrepancy may have implications for the MLE function worth exploring



## 1. Gompertz functions

library(flexsurv)
source("gompertz_functions.R") ## has functions like rgomp_mode()



## 2. Read in numident data (already truncated) from CenSoc demo file

library(data.table) ## to read in numident data
dt = fread("./censoc_numident_demo_v1/censoc_numident_demo_v1.csv")
dt1 <- dt[sex == "Male" & dyear %in% 1988:2005 & byear %in% 1918] ## restrict the sample
table(dt1$death_age, exclude = F) ## frequency table for estimating M
hist(dt1$death_age)  ## note that death_age values are floored integers
n = 1436 ## post-filter count of males
M.guess = 79 ## based on visual inspection of freq table and histogram
u = 88 ## upper bound, based on n = 1436  # upper bound needs to be one integer higher than observed max age at death
l = 69 ## lower bound, based on n = 1436
## beta is typically around 0.1
beta.guess = 0.1
y <- dt1[, death_age] ## here we use actual ages at death instead of simulated ages
yt = y ## actual ages are already truncated, so y is also yt



## 3. MLE estimation with known truncation

minusLogLik_1 =  function(yt, p, u, l)  ## p is a vector of M and log(beta)
{
    M = exp(p[1])
    beta = exp(p[2])
    n = length(yt)
    ## denom : F(u) - F(l)
    denom = pgomp_mode(u, M = M, beta = beta) - pgomp_mode(l, M = M, beta = beta)
    ## L = f1/denom * f2/denom ... = f1*f2*f3 / denom^n
    ## LL = sum(log(fi)) - n * log(denom)
    LL = sum(dgomp_mode(yt, M = M, beta = beta, log = TRUE)) - n* log(denom)
    mLL = -LL
    return(mLL)
}



## 4. Maximize the likelihood and find Hessian estimated standard errors

## starting values
p.start = c(log.M = log(M.guess),  
            log.beta = log(beta.guess))

## optimizer (to find the smallest minusLogLik using p.start parameters)
fit1 = optim(par = p.start, fn = minusLogLik_1, yt = yt, u = u, l = l)

## get out estimates
log.est = fit1$par
est = round(exp(log.est),3)
names(est) = c("M.hat", "beta.hat")
guess <- c("M.guess" = M.guess, "beta.guess" = beta.guess)
print(cbind(est, guess))

## Hessian estimated std err
set.seed(123)
fit1.hess = optim(par = p.start, fn = minusLogLik_1, yt = yt, u = u, l = l,
                  hessian = TRUE)
fit <- fit1.hess

## fisher_info<-solve(-fit$hessian)
H = fit$hessian
fisher_info = solve(H)
sigma.hat <- sqrt(diag(fisher_info))
upper <- fit$par + 1.96 * sigma.hat
lower <- fit$par - 1.96 * sigma.hat
interval <- data.frame(value = fit$par, upper = upper, lower = lower)
interval
exp(interval[1:2,])

