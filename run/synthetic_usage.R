set.seed(2)
ts <- c(1:20)
l <- 4
mvar <- 0.4
dist_m <- as.matrix(dist(ts,"manhattan",T,T))
cov_m <- apply(dist_m, c(1,2), function(d) mvar*exp(-0.5*(d/l)^2))

mean_func <- function(x) 0.2+(0.2/2)*max(min(x-8,2),0) + max(min(x-10, 4),0)*(-0.25/4) + max(min(x-14, 1),0)*(0.01/1) + max(min(x-16, 1),0)*(-0.05/1)
logi <- function(x) 1/(1+exp(-x))

mu <- sapply(ts, mean_func)
log_odds_mu <- sapply(mu, function(x) log(x/(1-x)))

phi <- rnorm(20, 0, 1)
M<-chol(cov_m)

logit_usage <- log_odds_mu + M %*% phi
usage_df <- data.frame(time=ts, usage=sapply(logit_usage,logi))
write.table(usage_df, "./../inst/extdata/synthetic_usage.csv", sep=",")