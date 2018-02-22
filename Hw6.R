#title : 405-Computational Methods HW6
#author: Ming-Hao Yu
#date  : 2018-02-21

FixStrikeLookbackOption <- function(s0, k, sigma, r, t, n, paths, seed, option="Call") {
    #s0: stock price, sigma: volatility, r: intereset rate, t: matrity time
    #n: steps, paths: times of simulations, seed: random generate seed, payoff: Call/Put
    set.seed(seed)
    dt <- t/n
    w1 <- sqrt(dt)*rnorm(paths/2*n)
    w <- c(w1, -w1)
    dwTable <- matrix(nrow=paths, ncol=n)
    stockPriceTable <- matrix(nrow=paths, ncol=n+1, s0)
    smax <- smin <- vector(length=paths)
    payoff <- vector(length=paths)
    
    #generate dw
    for(i in 1:paths) {
        dwTable[i, ] <- w[((i-1)*n+1):(i*n)]
    }
    #use simulated random walks to generate simulated stock prices
    for(i in 1:n) {
        stockPriceTable[, i+1] <- stockPriceTable[, i] + r*stockPriceTable[, i]*dt + sigma*stockPriceTable[, i]*dwTable[, i]
    }
    
    #get Smax at each path
    for(i in 1:paths) {
        smax[i] <- max(stockPriceTable[i,])
        smin[i] <- min(stockPriceTable[i,])
    }
    
    #calculate option prices
    if(toupper(option)=="CALL") {
        for(i in 1:paths) {
            payoff[i] <- max(smax[i]-k,0)
        }
    }
    else if(toupper(option)=="PUT") {
        for(i in 1:paths) {
            payoff[i] <- max(k-smin[i],0)
        }
    }
    else{cat("[No", option, "type of option.")}
    
    return(exp(-r*t)*mean(payoff))
}  

Problem1 <- function(seed) {
    s0 <- 98
    k <- 100
    t <- 1
    r <- 0.03
    paths <- 10000 #times of simulation
    n <- 200 #steps
    sigma <- seq(0.12,0.48, by=0.04)
    call <- vector(length=length(sigma))
    put <- vector(length=length(sigma))
    set.seed(seed)
    randseed <- round(runif(length(sigma))*10000)
    
    
    for(i in 1:length(sigma)) {
        call[i] <- FixStrikeLookbackOption(s0=s0, k=k, sigma=sigma[i], r=r, t=t, n=n, paths=paths, seed=randseed[i], option="Call")
        put[i] <- FixStrikeLookbackOption(s0=s0, k=k, sigma=sigma[i], r=r, t=t, n=n, paths=paths, seed=randseed[i], option="Put")
    }
    plot(y=call, x=sigma, main="Lookback Call Options", xlab="Volatility", ylab="Option Price", col="red", type="l")
    plot(y=put, x=sigma, main="Lookback Put Options", xlab="Volatility", ylab="Option Price", col="blue", type="l")
    return(list(Call=call, Put=put))
}

Problem2 <- function(lambda1, lambda2, T) {
    V0 <- 20000
    L0 <- 22000
    mu <- -0.1
    sigma <- 0.2
    gamma <- -0.4
    r0 <- 0.02
    delta <- 0.25
    alpha <- 0.7
    epsilon <- 0.95
    beta <- (epsilon-alpha)/T
    R <- r0+delta*lambda2
    r <- R/12
    n <- T*12
    PMT <- L0*r/(1-1/(1+r)^n)
    a <- PMT/r
    b <- PMT/(r*(1+r)^n)
    c <- 1+r
    
    set.seed(0)
    randomseed <- round(runif(3)*10000) #generate random seed from random variable distrbutions 
    paths <- 10000 #times of simulation
    steps <- 200 #steps
    dt <- T/steps
    
    set.seed(randomseed[1])
    w <- sqrt(dt)*rnorm(paths*steps)
    set.seed(randomseed[2])
    J <- rpois(n=paths*steps, lambda1*dt)
    set.seed(randomseed[3])
    N <- rpois(n=paths*steps, lambda2*dt)
    dwTable <- matrix(nrow=paths, ncol=steps)
    dJTable <- matrix(nrow=paths, ncol=steps)
    dNTable <- matrix(nrow=paths, ncol=steps)
    VTable <- matrix(nrow=paths, ncol=steps+1, V0)
    L <- rep(L0, steps+1)
    qt <- rep(epsilon,steps+1)
    
    Q <- vector(length=paths)
    S <- vector(length=paths)
    payoff <- vector(length=paths)
    exerciseTime <- vector(length=paths)
    #generate dw, dJ, dN table
    for(i in 1:paths) {
        dwTable[i, ] <- w[((i-1)*steps+1):(i*steps)]
        dJTable[i, ] <- J[((i-1)*steps+1):(i*steps)]
        dNTable[i, ] <- N[((i-1)*steps+1):(i*steps)]
    }
    
    #simulate each steps of V process: dVt=mu*Vt-*dt+sigma*Vt-*dWt+gamma*Vt-*dJt
    for(i in 1:steps) {
        VTable[, i+1] <- VTable[, i] + mu*VTable[, i]*dt + sigma*VTable[, i]*dwTable[, i] + gamma*VTable[, i]*dJTable[, i]
        L[i] <- a-b*c^(12*(i-1)*dt)
        qt[i] <- alpha+beta*(i-1)*dt
    }
    
    qt[1] <- alpha
    L[1] <- L0
    L[steps+1] <- 0
    qt[steps+1] <- epsilon
    
    #calculate Q and S
    for(i in 1:paths) {
        #calculate exercise time
        exerciseT1 <- which(VTable[i,] <= qt*L)
        exerciseT2 <- which(dNTable[i,]>0)
        Q[i] <- ifelse(length(exerciseT1)==0, steps, min(exerciseT1))
        S[i] <- ifelse(length(exerciseT2)==0, steps, min(exerciseT2))
        
        #calculate payoff
        exerciseTime[i] <- min(Q[i],S[i])
        if(exerciseTime[i]==200) {
            payoff[i] <- 0
        }
        else if(Q[i]<S[i]) {
            payoff[i] <- exp(-r0*Q[i]*dt)*max(L[Q[i]]-epsilon*VTable[i,Q[i]], 0)
        }
        else if(Q[i]==S[i]){
            payoff[i] <- max(max(L[Q[i]]-epsilon*VTable[i,Q[i]], 0), abs(L[S[i]]-epsilon*VTable[i,S[i]]))
        }
        else {
            payoff[i] <- exp(-r0*S[i]*dt)*abs(L[S[i]]-epsilon*VTable[i,S[i]])
        }
    }
    
    D <- mean(payoff)
    Prob <- length(which(exerciseTime!=200))/paths
    Et <- mean(exerciseTime[which(exerciseTime!=200)]/steps*T)
    return(list(D=D, Prob=Prob, Et=Et))
}

Problem1(0)
Problem2(lambda1=0.2, lambda2=0.4, T=5)

lambda1 <- seq(0.05,0.4, by=0.05)
lambda2 <- seq(0, 0.8, by=0.1)
T <- seq(3,8,by=1)

price1 <- matrix(nrow=length(lambda2), ncol=length(T))
price2 <- matrix(nrow=length(lambda1), ncol=length(T))
Prob1 <- matrix(nrow=length(lambda2), ncol=length(T))
Prob2 <- matrix(nrow=length(lambda1), ncol=length(T))
Et1 <- matrix(nrow=length(lambda2), ncol=length(T))
Et2 <- matrix(nrow=length(lambda1), ncol=length(T))
colnames(price1) <- T
rownames(price1) <- lambda2
colnames(price2) <- T
rownames(price2) <- lambda1
colnames(Prob1) <- T
rownames(Prob1) <- lambda2
colnames(Prob2) <- T
rownames(Prob2) <- lambda1
colnames(Et1) <- T
rownames(Et1) <- lambda2
colnames(Et2) <- T
rownames(Et2) <- lambda1
for(i in 1:length(lambda2)) {
    for(j in 1:length(T)) {
        ans <- Problem2(lambda1=0.2, lambda2=lambda2[i], T[j])
        price1[i, j] <- ans$D
        Prob1[i, j] <- ans$Prob
        Et1[i, j] <- ans$Et
    }
}
for(i in 1:length(lambda1)) {
    for(j in 1:length(T)) {
        ans <- Problem2(lambda1=lambda1[i], lambda2=0.4, T[j])
        price2[i, j] <- ans$D
        Prob2[i, j] <- ans$Prob
        Et2[i, j] <- ans$Et
    }
}

###### Plots #####
### Price vs T
ylim <- c(min(price1), max(price1))
for(i in 1:length(lambda2)) {
    plot(y=price1[i,], x=T, main="Option price, lambda1=0.4", xlab="Maturity T", ylab="Price", col=i, type="l", ylim=ylim)
    par(new=TRUE)
}
par(new=FALSE)
legend("bottomright", legend=lambda2, title="lambda2", col=seq(1,length(lambda2)), lwd=2.5, cex=0.8)

ylim <- c(min(price2), max(price2))
for(i in 1:length(lambda1)) {
    plot(y=price2[i,], x=T, main="Option price, lambda2=0.2", xlab="Maturity T", ylab="Price", col=i, type="l", ylim=ylim)
    par(new=TRUE)
}
par(new=FALSE)
legend("bottomright", legend=lambda1, title="lambda1", col=seq(1,length(lambda1)), lwd=2.5, cex=0.8)

### Probability vs T
ylim <- c(min(Prob1), max(Prob1))
for(i in 1:length(lambda2)) {
  plot(y=Prob1[i,], x=T, main="Default probability, lambda1=0.4", xlab="Maturity T", ylab="Probability", col=i, type="l", ylim=ylim)
  par(new=TRUE)
}
par(new=FALSE)
legend("bottomright", legend=lambda2, title="lambda2", col=seq(1,length(lambda2)), lwd=2.5, cex=0.8)

ylim <- c(min(Prob2), max(Prob2))
for(i in 1:length(lambda1)) {
  plot(y=Prob2[i,], x=T, main="Default probability, lambda2=0.2", xlab="Maturity T", ylab="Probability", col=i, type="l", ylim=ylim)
  par(new=TRUE)
}
par(new=FALSE)
legend("bottomright", legend=lambda1, title="lambda1", col=seq(1,length(lambda1)), lwd=2.5, cex=0.8)

###  Estimated Exercise Time vs T
ylim <- c(min(Et1), max(Et1))
for(i in 1:length(lambda2)) {
  plot(y=Et1[i,], x=T, main="Estimated Exercise Time, lambda1=0.4", xlab="Maturity T", ylab="Year", col=i, type="l", ylim=ylim)
  par(new=TRUE)
}
par(new=FALSE)
legend("bottomright", legend=lambda2, title="lambda2", col=seq(1,length(lambda2)), lwd=2.5, cex=0.8)

ylim <- c(min(Et2), max(Et2))
for(i in 1:length(lambda1)) {
  plot(y=Et2[i,], x=T, main="Estimated Exercise Time, lambda2=0.2", xlab="Maturity T", ylab="Year", col=i, type="l", ylim=ylim)
  par(new=TRUE)
}
par(new=FALSE)
legend("bottomright", legend=lambda1, title="lambda1", col=seq(1,length(lambda1)), lwd=2.5, cex=0.8)
