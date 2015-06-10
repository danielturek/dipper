

#################################################################
###########  llFunction for mark-recapture histories:  ##########
###########     phi: constant / time dependence        ##########
###########     p:   constant / time dependence        ##########
#################################################################


llFunction_MR_CJS <- nimbleFunction(
    
    setup = function(model, dataNode, phiNode, pNode) {
        y <- as.matrix(model[[dataNode]])
        nind  <- dim(y)[1]
        k     <- dim(y)[2]
        first <- apply(y, 1, function(hist) min(which(hist == 1)))
        last  <- apply(y, 1, function(hist) max(which(hist == 1)))
        
        phiNode <- model$expandNodeNames(phiNode)
        pNode   <- model$expandNodeNames(pNode)
        needToRepPhi <- if(length(phiNode) == 1) TRUE else { if(length(phiNode) != k) stop('phiNode has wrong length'); FALSE }
        needToRepP   <- if(length(pNode)   == 1) TRUE else { if(length(pNode)   != k) stop('pNode has wrong length');   FALSE }
        phiVec <- rep(0, k)
        pVec   <- rep(0, k)
    },
    
    run = function() {
        declare(phiVec, double(1, k))
        declare(pVec,   double(1, k))
        if(needToRepPhi) { phiValue <- model[[phiNode]];   for(i in 1:k) phiVec[i] <<- phiValue }
        else             { getValues(phiVec, model, phiNode) }
        if(needToRepP)   { pValue   <- model[[pNode]];     for(i in 1:k) pVec[i]   <<- pValue   }
        else             { getValues(pVec,   model, pNode)   }
        
        declare(chi, double(1, k))
        chi[k] <- 1
        for(i in 1:(k-1))     { chi[k-i] <- 1 - phiVec[k-i] + phiVec[k-i]*(1-pVec[k-i+1])*chi[k-i+1] }
        logChi <- log(chi)
        
        declare(llVec, double(1, nind))
        
        for(ind in 1:nind) {
            ll <- 0
            if(first[ind] < last[ind]) {
                for(t in (first[ind]+1) : last[ind] )     { ll <- ll + log(phiVec[t-1]) + y[ind,t]*log(pVec[t]) + (1-y[ind,t])*log(1-pVec[t]) } }
            ll <- ll + logChi[last[ind]]
            llVec[ind] <- ll
        }
        logLikelihood <- sum(llVec)
        returnType(double())
        return(logLikelihood)
    }
)




#################################################################
###########  llFunction for mark-recapture histories:  ##########
###########     phi: constant                          ##########
###########     p:   constant                          ##########
#################################################################

llFunction_MR_phiConst_pConst <- nimbleFunction(

    setup = function(model, dataNode, phiNode, pNode) {
        y <- as.matrix(model[[dataNode]])    ## as.matrix() necessary, to get dimensions correct
        y <- 2 - y    ## changes y=1 (alive) to state=1, and y=0 (dead) to state=2
        nind  <- dim(y)[1]   ## number of individual MR histories
        k     <- dim(y)[2]   ## number of sighting occasions
        first <- apply(y, 1, function(MRhistory) min(which(MRhistory == 1)))
        nStates <- 2
        piXfirst <- matrix(c(1, 0))   ## prior distribution of state variable at first encounter
    },
    
    run = function() {
        declare(first, double(1, nind))    ## delcare() is necessary here, to ensure 'first' has dimension = 1
        
        if(model[[phiNode]] < 0) return(NaN)   ## return(-Inf) instead ????
        if(model[[phiNode]] > 1) return(NaN)   ##
        if(model[[pNode]]   < 0) return(NaN)   ##
        if(model[[pNode]]   > 1) return(NaN)   ##
        
        declare(Tmat, double(2, c(nStates, nStates)))   ## state transition matrix
        Tmat[1, 1] <- model[[phiNode]]
        Tmat[2, 1] <- 1 - model[[phiNode]]
        Tmat[1, 2] <- 0
        Tmat[2, 2] <- 1
        
        declare(Zmat, double(2, c(nStates, nStates)))   ## observation process matrix
        Zmat[1,1] <- model[[pNode]]
        Zmat[2,1] <- 1 - model[[pNode]]
        Zmat[1,2] <- 0
        Zmat[2,2] <- 1
        
        declare(Lvec, double(1, nind))   ## vector of conditional likelihood values
        
        for(ind in 1:nind) {
            piX <- piXfirst
            for(t in first[ind]:k) {
                declare(Zslice, double(2, c(2,1)))
                for(iState in 1:nStates)     { Zslice[iState, 1] <- Zmat[y[ind, t], iState] }    ## 'slice' of Z matrix relevant to y[ind, t]
                ## initialize or update cumulative likelihood
                if(t == first[ind])     { Lvec[ind] <- 1    ## first encounter (conditioned on y=1), initialize L <- 1
                } else                  { Lvec[ind] <- Lvec[ind] * sum(piX * Zslice) }
                piXstar <- piX * Zslice / sum(piX * Zslice)   ## update distribution of state vector, conditional on observation y[ind, t]
                piX <- Tmat %*% piXstar   ## propagate state distribution to time (t+1)
            }
            ##print('Lvec[ind] = ', Lvec[ind])
        }
        L <- prod(Lvec)
        returnType(double())
        return(log(L))
    }
    #where = getNamespace('nimble')
)


#################################################################
###########  llFunction for mark-recapture histories:  ##########
###########     phi: constant                          ##########
###########     p:   individual heterogeneity          ##########
#################################################################


llFunction_MR_phiConst_pHeter <- nimbleFunction(

    setup = function(model, dataNode, phiNode, pNodeVector, capture_hist_indices = NULL) {
        y <- as.matrix(model[[dataNode]])    ## as.matrix() necessary, to get dimensions correct
        y <- 2 - y    ## changes y=1 (alive) to state=1, and y=0 (dead) to state=2
        if(!is.null(capture_hist_indices)) {
            y <- y[capture_hist_indices, , drop = FALSE]
            pNodeVector <- model$expandNodeNames(pNodeVector)[capture_hist_indices]
        }
        nind  <- dim(y)[1]   ## number of individual MR histories
        k     <- dim(y)[2]   ## number of sighting occasions
        first <- apply(y, 1, function(MRhistory) min(which(MRhistory == 1)))
        nStates <- 2
        piXfirst <- matrix(c(1, 0))   ## prior distribution of state variable at first encounter
    },
    
    run = function() {
        declare(first, double(1, nind))    ## delcare() is necessary here, to ensure 'first' has dimension = 1
        
        if(model[[phiNode]] < 0) return(NaN)   ## return(-Inf) instead ????
        if(model[[phiNode]] > 1) return(NaN)   ##
        
        declare(Tmat, double(2, c(nStates, nStates)))   ## state transition matrix
        Tmat[1, 1] <- model[[phiNode]]
        Tmat[2, 1] <- 1 - model[[phiNode]]
        Tmat[1, 2] <- 0
        Tmat[2, 2] <- 1
        
        declare(pValues, double(1, nind))
        getValues(pValues, model, pNodeVector)
        
        declare(Lvec, double(1, nind))   ## vector of conditional likelihood values
        
        for(ind in 1:nind) {
            if(pValues[ind] < 0) return(NaN)   ## return(-Inf) instead ????
            if(pValues[ind] > 1) return(NaN)   ##
            
            declare(Zmat, double(2, c(nStates, nStates)))   ## observation process matrix
            Zmat[1,1] <- pValues[ind]
            Zmat[2,1] <- 1 - pValues[ind]
            Zmat[1,2] <- 0
            Zmat[2,2] <- 1
            
            piX <- piXfirst
            
            for(t in first[ind]:k) {
                declare(Zslice, double(2, c(2,1)))
                for(iState in 1:nStates)     { Zslice[iState, 1] <- Zmat[y[ind, t], iState] }    ## 'slice' of Z matrix relevant to y[ind, t]
                ## initialize or update cumulative likelihood
                if(t == first[ind])     { Lvec[ind] <- 1    ## first encounter (conditioned on y=1), initialize L <- 1
                } else                  { Lvec[ind] <- Lvec[ind] * sum(piX * Zslice) }
                piXstar <- piX * Zslice / sum(piX * Zslice)   ## update distribution of state vector, conditional on observation y[ind, t]
                piX <- Tmat %*% piXstar   ## propagate state distribution to time (t+1)
            }
        }
        L <- prod(Lvec)
        returnType(double())
        return(log(L))
    }
    #where = getNamespace('nimble')
)



