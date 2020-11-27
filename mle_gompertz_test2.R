## mle_gompertz_test2.R

## same as previous test, but now we account for covariates

## Note: we use proportional hazards Gompertz model
##    h(x)[i] = alpha * exp(beta * x) * exp(b * x[i])
##    h(y)[i] = alpha * exp(beta * y) * exp(b * x[i])
## where exp(b * x[i]) is the proportional effect of x on hazards



## 1. Gompertz functions

library(flexsurv)
source("gompertz_functions.R")



## 2. Read in numident data (already truncated) from CenSoc demo file and clean up incwage covariate

library(data.table) ## to read in numident data
dt = fread("./censoc_numident_demo_v1/censoc_numident_demo_v1.csv")
dt[, inc := incwage]
dt[incwage == 0, inc := NA] ## 0 doesn't mean no earnings, it means missing for various reasons
dt[incwage == 999999, inc := NA] ## 999999 means not eligible for question for various reasons
dt$log_inc <- log(dt$inc) ## take the log of "clean" inc

dt2 <- dt[sex == "Male" & dyear %in% 1988:2005 & byear %in% 1918 & !is.na(log_inc)] ## create sample
table(dt2$death_age, exclude = F) ## frequency table for estimating M
hist(dt2$death_age) ## note that death_age values are floored integers
hist(dt$log_inc) ## distribution skews left

n = 990 ## post-filter count of males
M.guess = 79 ## based on visual inspection of freq table and histogram
u = 88 ## upper bound, based on n = 990; upper bound needs to be one integer higher than observed max age at death
l = 69 ## lower bound, based on n = 990
beta.guess = 1/10

x = dt2$log_inc ## covariate is the log of wage income
#x = scale(dt2$log_inc) ## covariate is the log of wage income after scaling

b.guess = .6 ## the coefficient (effect) we will try to estimate
alpha = 10^-4.5 ## alpha tends to be between 10^-6 and 10^-3
rate.vec = alpha * exp(x * b.guess) ## proportional hazard parameters for our observed sample

M.vec <- getMode(alpha = rate.vec, beta = beta.guess) ## in terms of mode
#M.vec <- dt2[, death_age]

## we use actual ages at death instead of simulated ages
y <- dt2[, death_age] 
## truncated death ages
yt = y ## actual ages are already truncated, so y is also yt
## corresponding truncated covariate
xt = x



## 3. MLE estimation with known truncation

minusLogLik_2 =  function(yt, xt, p, u, l) ## p is a vector of M and log(beta)
{
    M = exp(p[1])
    beta = exp(p[2])
    b = p[3] ## no exp()
    alpha = getAlpha(M = M, beta = beta) ## may be thought of as alpha.hat and beta.hat
    rate.vec = alpha * exp(xt * b)  ## may be thought of as rate.vec.hat
    M.vec <- getMode(alpha = rate.vec, beta = beta) 
    ## denom : F(u) - F(l)
    denom = pgomp_mode(u, M = M.vec, beta = beta) - pgomp_mode(l, M = M.vec, beta = beta)
    ## replace 0s with a very small number for numerical stability (before taking logs)
    eps = 10^-6
    denom[denom == 0] <- eps
    ## L = f1/denom * f2/denom ... = f1*f2*f3 / denom^n
    ## LL = sum(log(fi)) - sum( log(denom))
    LL = sum(dgomp_mode(yt, M = M.vec, beta = beta, log = TRUE)) - sum(log(denom))
    mLL = -LL
    return(mLL)
}



## 4. Maximize the likelihood and find Hessian estimated standard errors

## starting values
p.start = c(log.M = log(M.guess * .9), ## kind of close to guessed values
            log.beta = log(beta.guess * 1.2),
            b = b.guess * .5)

## optimizer
fit2 = optim(par = p.start, fn = minusLogLik_2, yt = yt, xt = xt, u = u, l = l)

## get out estimates
log.est = fit2$par
est = round(c(exp(log.est[1:2]), log.est[3]),3)
names(est) = c("M.hat", "beta.hat", "b.hat")
guess <- c("M.guess" = M.guess, "beta.guess" = beta.guess, b.guess = b.guess)
print(cbind(est, guess))

## Hessian estimated std err
set.seed(123)
fit2.hess = optim(par = p.start, fn = minusLogLik_2, yt = yt, xt = xt, u = u, l = l,
                  hessian = TRUE)
fit <- fit2.hess

## fisher_info<-solve(-fit$hessian)
H = fit$hessian
fisher_info = solve(H)
sigma.hat <- sqrt(diag(fisher_info))
upper <- fit$par + 1.96 * sigma.hat
lower <- fit$par - 1.96 * sigma.hat
interval <- data.frame(value = fit$par, upper = upper, lower = lower)
interval
exp(interval[1:2,])

