


load('dipperData.RData')



## optionally, truncate data
ind <- 1:3;   nind<-length(ind);   first<-first[ind];   y<-y[ind,,drop=FALSE];   x_init<-x_init[ind,,drop=FALSE]



code_dipper <- nimbleCode({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for (i in 1:nind) {
        x[i, first[i]] <- 1
        for (t in (first[i] + 1):k) {
            mu_x[i, t] <- phi * x[i, t-1]
            mu_y[i, t] <- p * x[i, t]
            x[i, t] ~ dbin(mu_x[i, t], 1)
            y[i, t] ~ dbin(mu_y[i, t], 1)
        }
    }
})
constants <- list(k=k, nind=nind, first=first)
data      <- list(y=y)
inits     <- list(phi=0.6, p=0.9, x=x_init)



code_dipperPH <- nimbleCode({
    phi ~ dunif(0, 1)
    pmean ~ dunif(0, 1)
    plogitmean <- log(pmean/(1 - pmean))
    psigma ~ dunif(0, 5)
    ptau <- pow(psigma, -2)
    for (i in 1:nind) {
        epsilon[i] ~ dnorm(0, ptau)
        logit(p[i]) <- plogitmean + epsilon[i]
        x[i, first[i]] <- 1
        for (t in (first[i] + 1):k) {
            mu_x[i, t] <- phi * x[i, t - 1]
            mu_y[i, t] <- p[i] * x[i, t]
            x[i, t] ~ dbin(mu_x[i, t], 1)
            y[i, t] ~ dbin(mu_y[i, t], 1)
        }
    }
})
constants <- list(k=k, nind=nind, first=first)
data      <- list(y=y)
inits     <- list(phi=0.6, pmean=0.9, psigma=1, x=x_init)   ## WinBUGS sensitive to psigma initial value



code_dipperSeasonal <- nimbleCode({
    phi_flood ~ dunif(0,1)
    phi_non   ~ dunif(0,1)
    p         ~ dunif(0,1)
    phi[1] <- phi_non
    phi[2] <- phi_flood
    phi[3] <- phi_flood
    phi[4] <- phi_non
    phi[5] <- phi_non
    phi[6] <- phi_non
    phi[7] <- 1
    for(i in 1:nind) {
        x[i, first[i]] <- 1
        for(t in (first[i]+1):k) {
            mu_x[i,t] <- phi[t-1] * x[i,t-1]
            mu_y[i,t] <- p        * x[i,t]
            x[i,t] ~ dbern(mu_x[i,t])
            y[i,t] ~ dbern(mu_y[i,t])
        }
    }
})
constants <- list(k=k, nind=nind, first=first)
data      <- list(y=y)
inits     <- list(phi_flood=0.6, phi_non=0.6, p=0.9, x=x_init)





