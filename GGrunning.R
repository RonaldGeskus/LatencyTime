
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
