## mle_gompertz_simu.R

## josh g. (Nov 9, 2020)

## Overview: we demonstrate that the truncated gompertz is possible to fit successfully to
## simulated gompertz data

## As extra-credit, we do with a covariate

## And as extra-extra credit, we estimate SE of our estimates using Hessian

## Notes:
## * We use mode M parameterization of Gompertz
## * Not sure which optimizer works best but I just use default
## * Problems with optimization often are coding problems with the objective function





## for birth cohort of 1918, OLS is estimating: beta = d(Y | Y %in% 70:87) / dX
## Achim-method : beta_a = d(Y | Y <= 87) / dX

## Could we invert the problem? Not completely, because we're dealing with survival ... and that affects the way we think about it

## Achim-method-inverted? : beta_a = d(Y | Y >= 70) / dX
## What do we really want? Something like beta_true = d(Y | Y >= 65) / dX

## For example, in parametric MLE:
## L = f(y | alpha(X), beta(X)) / [F(u | alpha(X), beta(X)) - F(l | alpha(X), beta(X))]
## then with alpha.hat and beta.hat and X, I would estimate d(Y(X) | Y >= 65) / dX







## Contents:
## 1. Gompertz functions
## 2. Simulate data (without covariate) with truncation
## 3. ML estimates
## 4. Simulate data with a covariate (with truncation) and re-estimate
## 5. Hessian example for estimated standard errors

setwd("~/Google Drive/berkeley/CenSoc/gompertz")
library(txtplot) ## for ascii plots within terminal


## 1. Gompertz functions

library(flexsurv)
source("gompertz_functions.R") ## has functions like rgomp_mode()


## 2. Simulate data (without covariate) with truncation

set.seed(314)
n = 1000
M.true = 88
beta.true = 1/10
u = 90 ## upper bound
l = 80 ## lower bound
y = rgomp_mode(n = n, M = M.true, beta = beta.true) ## untruncated ages at death
yt = y[l < y & y < u] ## truncated ages at death


## 3. MLE estimation with known truncation

minusLogLik_1 =  function(yt, p, u, l)
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

## Maximize the likelihood (by minimizing the minusLogLik)

## starting values
p.start = c(log.M = log(M.true * .8), ## kind of close to true values
            log.beta = log(beta.true * 1.2))

## optimizer
fit1 = optim(par = p.start, fn = minusLogLik_1, yt = yt, u = u, l = l)

## get out estimates
log.est = fit1$par
est = round(exp(log.est),3)
names(est) = c("M.hat", "beta.hat")
true <- c("M.true" = M.true, "beta.true" = beta.true)
print(cbind(est, true))
##             est true
## M.hat    87.313 88.0
## beta.hat  0.088  0.1

## pretty good match (if you increase sample size to 10000, could probably do even better) or could simulate this smaller sample many times

## Note: you can use Hessian form optimizer to get estimated variance (and SE) of estimates



## 4. Simulate data with a covariate (with truncation) and re-estimate

## Note: we use proportional hazards Gompertz model
## h(x)[i] = alpha * exp(beta * x) * exp(b * x[i])
## h(y)[i] = alpha * exp(beta * y) * exp(b * x[i])
## where exp(b * x[i]) is the proportional effect of x on hazards

## Simulate with n obs
## set.seed(23)
n = 10000 ## bigger sample
x = rnorm(n) ## just let the covariate be normal(0,1)
b.true = .6 ## the coefficient (effect) we will try to estimate
alpha = 10^-4.5
rate.vec = alpha * exp(x * b.true) ## proportional hazard parameters for our simulated sample
M.vec <- getMode(alpha = rate.vec, beta = beta.true) ## in terms of mode
## now simulate, taking advantage of ability to use vectors as arguments for parameter values
y = rgomp_mode(n = n, M = M.vec, beta = beta.true)
## truncated death ages
yt = y[l < y & y < u] ## truncated ages at death
## corresponding truncated covariate
xt = x[l < y & y < u]


minusLogLik_2 =  function(yt, xt, p, u, l)
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

## starting values
p.start = c(log.M = log(M.true * .9), ## kind of close to true values
            log.beta = log(beta.true * 1.2),
            b = b.true * .5)

## optimizer
fit2 = optim(par = p.start, fn = minusLogLik_2, yt = yt, xt = xt, u = u, l = l)

## get out estimates
log.est = fit2$par
est = round(c(exp(log.est[1:2]), log.est[3]),3)
names(est) = c("M.hat", "beta.hat", "b.hat")
true <- c("M.true" = M.true, "beta.true" = beta.true, b.true = b.true)
print(cbind(est, true))
##             est true
## M.hat    78.982 88.0
## beta.hat  0.096  0.1
## b.hat     0.606  0.6
## very similar values


## 5. Hessian example for estimated standard errors

set.seed(123)
fit2.hess = optim(par = p.start, fn = minusLogLik_2, yt = yt, xt = xt, u = u, l = l,
                  hessian = TRUE)
##https://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
fit <- fit2.hess

## fisher_info<-solve(-fit$hessian)
H = fit$hessian
fisher_info = solve(H)
sigma.hat <- sqrt(diag(fisher_info))
upper <- fit$par + 1.96 * sigma.hat
lower <- fit$par - 1.96 * sigma.hat
interval <- data.frame(value = fit$par, upper = upper, lower = lower)
interval
##               value      upper      lower
## log.M     4.3798717  4.4337155  4.3260278
## log.beta -2.3445257 -1.8751319 -2.8139195
## b         0.6063259  0.8038096  0.4088421

exp(interval[1:2,])
##                value      upper       lower
## log.M    79.82778855 84.2438452 75.64322128
## log.beta  0.09589267  0.1533347  0.05996948

## ok, this all looks reasonable but I'm not 100% sure
## and I doubt that exp(interval) is really exactly correct
