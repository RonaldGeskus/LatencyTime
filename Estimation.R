#####################################################################################
### Estimation of latency time (data with strict exposure windows)
### Code for generalized gamma distribution
#####################################################################################
library(rjags)
library(ggplot2)
library(tidyr)
library(coda)

## data preparation
latency <- read.csv("Data/IntervalsStrict.csv")
## Subset of exposure windows < 5 days
## latency_narrow <- latency |> subset(L1 - L0 < 5)

## prepare data for JAGS program
N <- length(latency$L0)
dat <- with(latency, list(ER = L1-L0, SL = R0-L0, SR = R1-L0, Trunc = Trunc-L0,
                          type_S = rep(1L,N), r_est = 0.1058827, r_se = 0.004466976, N = N,
                          Y = rep(NA,N)))
#####################################################################################
### exponential growth for infection time
#####################################################################################

## function to generate initial values
jags_inits_growth <- function() {
    E_init <- Y_init <- rep(NA, N)
    for (i in 1:N) {
        if (dat$SL[i] == 0 & dat$ER[i] == dat$SR[i]) {
            tmp_inits <- runif(2, min = 0, max = dat$ER[i])
            E_init[i] <- min(tmp_inits)
            Y_init[i] <- max(tmp_inits) - E_init[i]
        }
    else {
        E_init[i] <- runif(1, min = 0, max = min(dat$ER[i], dat$SL[i]))
        Y_init[i] <- runif(1, min = max(dat$ER[i], dat$SL[i]), max = dat$SR[i]) - E_init[i]
        }
    }
    E_star_init <- dat$ER - E_init
    rc <- list(E_star = E_star_init,  Y = Y_init)
    return(rc)
}

## run JAGS
res.genG <- jags.model("JAGS/GGammaGrowth.jag", data = dat, inits=jags_inits_growth,
                     n.adapt = 2000, quiet = FALSE, n.chains = 3)
update(res.genG, n.iter=8000)
res.genG <- jags.samples(res.genG, c("r", "lambda", "b", "mean", "perc50", "perc90", "perc95", "perc97.5", "perc99"), n.iter = 500000,  quiet = TRUE, thin=10)

## check convergence and summarize results
gg.running(res.genG, quantile=FALSE)
gg.running(res.genG, thin=10)

lapply(res.genG, FUN=function(x) summary(as.mcmc.list(x)))

#####################################################################################
### uniform distribution for infection time
#####################################################################################

## function to generate initial values
jags_inits_unif <- function() {
    E_init <- Y_init <- rep(NA, N)
    for (i in 1:N) {
        if (dat$SL[i] == 0 & dat$ER[i] == dat$SR[i]) {
            tmp_inits <- runif(2, min = 0, max = dat$ER[i])
            E_init[i] <- min(tmp_inits)
            Y_init[i] <- max(tmp_inits) - E_init[i]
        }
    else {
        E_init[i] <- runif(1, min = 0, max = min(dat$ER[i], dat$SL[i]))
        Y_init[i] <- runif(1, min = max(dat$ER[i], dat$SL[i]), max = dat$SR[i]) - E_init[i]
        }
    }
    rc <- list(E = E_init, b = runif(1, 0, 10),
               r = runif(1, 0, 10), lambda = runif(1, 0, 10), Y = Y_init)
    return(rc)
}

## run JAGS
res.genG <- jags.model("JAGS/GGammaUnif.jag", data = dat, inits=jags_inits_unif,
                     n.adapt = 2000, quiet = FALSE, n.chains = 3)
update(res.genG, n.iter=8000)
res.genG <- jags.samples(res.genG, c("r", "lambda", "b", "mean", "perc50", "perc90", "perc95", "perc97.5", "perc99"), n.iter = 500000,  quiet = TRUE, thin=10)

## check convergence and summarize results
gg.running(res.genG, quantile=FALSE)
gg.running(res.genG, thin=10)

lapply(res.genG, FUN=function(x) summary(as.mcmc.list(x)))

#####################################################################################
### uniform distribution for infection time
#####################################################################################

Est <- Estimate_doublIn(latency, infection_risk_distribution = "exp_growth",
                        exp_growth_rate = 0.1058827,  exp_growth_rate_SE = 0.004466976,
                        right_truncation = TRUE, iters = 500000, burnin_period = 10000,
                        thin = 10)





gg.running <- function(mcmc,  sels=NULL, method="JAGS", probs=c(0.025, 0.5, 0.975), thin=1,
                       quantile=TRUE){
    ## plots chains (if quantile=FALSE) or running quantiles (if quantile=TRUE) in ggplot2
    ## mcmc: list of mcarray's (JAGS)
    ## sels: JAGS: named list to select elements of the components in mcmc;

    ## chains: list, each chain one component, each component and array with n.iter rows and n.col parameters
    ## vars: the names of the parameters
    require(rjags)
    require(tidyr)
    theme_set(theme_bw())
    cquantile <- function(z, probs) {
        cquant <- matrix(0, nrow = length(z), length(probs))
        for (i in seq(along = z))
            cquant[i, ] <- quantile(z[1:i], probs = probs, names = FALSE)
        cquant <- as.data.frame(cquant)
        names(cquant) <- c("lower","median","upper")
        return(cquant)
    }
    Cumul <- vector(mode="list",length=length(names(mcmc)))
    names(Cumul) <- names(mcmc)
    n.chain <- dim(mcmc[[1]])["chain"]
    for(var in names(mcmc)){
        Cumul[[var]] <- as.mcmc.list(mcmc[[var]])
        L <- length(Cumul[[var]][[1]][,1])
        if(!is.null(sels)){
            var.names <- dimnames(Cumul[[var]][[1]])[[2]][sels[[var]]]
            for(i in 1:3){
                Cumul[[var]][[i]] <- Cumul[[var]][[i]][seq(1,L,by=thin),sels[[var]]]
            }
        } else {
            var.names <- dimnames(Cumul[[var]][[1]])[[2]]
            for(i in 1:n.chain){
                Cumul[[var]][[i]] <- Cumul[[var]][[i]][seq(1,L,by=thin),]
            }
        }
        if(is.null(dim(Cumul[[var]][[1]]))) nval <- length(Cumul[[var]][[1]]) else
                                            nval <- prod(dim(Cumul[[var]][[1]]))
        Cumul[[var]] <- data.frame(iter=seq(1,L,by=thin),
                    values=unlist(Cumul[[var]]),
                    chain=as.factor(rep(1:length(Cumul[[var]]), rep(nval,length(Cumul[[var]])))),
                    vars=rep(var.names, rep(length(seq(1,L,by=thin)),length(var.names))))
    }
    Cumul <- as.data.frame(do.call("rbind",Cumul))
    Cumul$vars <- factor(Cumul$vars,levels=unique(Cumul$vars))
    if(quantile){
        Cumul <- Cumul %>% group_by(vars,chain) %>% mutate(cquantile(values,c(0.025, 0.5, 0.975)))
        ggplot(Cumul) + geom_line(aes(iter,median,group=chain,color=chain))+
            geom_line(aes(iter,lower,group=chain,color=chain))+ geom_line(aes(iter,upper,group=chain,color=chain)) +
            facet_wrap(~vars,scales="free_y")
    } else {
        ggplot(Cumul) + geom_line(aes(iter,values,group=chain,color=chain),alpha=0.5) +
            facet_wrap(~vars,scales="free_y")+ theme(legend.position="null")
    }
}
