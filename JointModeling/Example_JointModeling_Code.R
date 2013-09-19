
#
# Code to fit a joint model for survival and longitudinal data
# using the MCECM algorithm. Note that due to MC variability in 
# the E-step the parameter values will usually just hover around
# the MLE. For strategies to obtain consistent estimates from
# these paths see Baines, Xu and Wang (2013) [not implemented here].
#

# Compile or not?
do_compile <- TRUE

library(MASS)
library(nlme)
library(survival)

###
#
# Functions:
#
# JMiterECM1 -- (Internal only) Runs a single JM ECM iteration
# JMinitECM  -- (Internal only) Computes "good" initial values for the algorithm
# JointModelECM -- (User) Fit joint model using ECM
#
###

#=============== ECM iteration ===============#
"JMiterECM1" <- function(theta_old,ModelInfo,MC,seed){

        # Get Old Estimates #
        beta_old = theta_old$beta
        Ysigma_old = theta_old$Ysigma
        BSigma_old = theta_old$BSigma
        gamma_old = theta_old$gamma
        alpha_old = theta_old$alpha
        lamb_old = theta_old$lamb

        # Get Data Information in the Model #
        ID = ModelInfo$ID # vector of length N #
        d = ModelInfo$d # vector of length n #
        Index = ModelInfo$Index # vector of length M #
        Index0 = ModelInfo$Index0 # vector of length n #
        Index1 = ModelInfo$Index1 # vector of length M #
        Index2 = ModelInfo$Index2 # vector of length n_u #
        U = ModelInfo$U # vector of length n_u #
        Y = ModelInfo$Y # vector of length N #
        X = as.matrix(ModelInfo$X) # N*ncx matrix #
        Z = as.matrix(ModelInfo$Z) # N*ncz matrix #
        Xtime = as.matrix(ModelInfo$Xtime) # n*ncx matrix #
        Xtime2 = as.matrix(ModelInfo$Xtime2) # M*ncx matrix #
        Ztime = as.matrix(ModelInfo$Ztime) # n*ncz matrix #
        Ztime2 = as.matrix(ModelInfo$Ztime2) # M*ncz matrix #
        Wtime = as.matrix(ModelInfo$Wtime) # n*ncw matrix #

        M = length(Index)
        N = length(Y)
        n_u = length(U)
        ncx = ncol(X)
        ncz = ncol(Z)
        ncw = ncol(Wtime)

        if(!is.null(seed)){
          set.seed(seed)
        }
        bi = mvrnorm(MC,rep(0,ncz),BSigma_old) # MC*ncz matrix #

        Ymu = as.vector(X%*%beta_old)+Z%*%t(bi) # N*MC matrix #
        logNorm = dnorm(Y,mean=Ymu,sd=sqrt(Ysigma_old),log=TRUE) # N*MC matrix #
        part1 = rowsum(logNorm,ID) # n*MC matrix #

        eta1 = as.vector(Xtime%*%beta_old)+Ztime%*%t(bi) # n*MC matrix #
        eta2 = as.vector(Xtime2%*%beta_old)+Ztime2%*%t(bi) # M*MC matrix #
        log.lamb = log(lamb_old[Index0])
        log.lamb[is.na(log.lamb)] = 0
        log.hazard = log.lamb+as.vector(Wtime%*%gamma_old)+alpha_old*eta1 # n*MC matrix #
        eta.s = as.vector(Wtime%*%gamma_old)[Index]+alpha_old*eta2 # M*MC matrix #
        log.survival = -rowsum(lamb_old[Index1]*exp(eta.s),Index) # n*MC matrix #
        part2 = d*log.hazard+log.survival # n*MC matrix #

        total = exp(part1+part2) # n*MC matrix #
        deno = as.vector(rowMeans(total)) # vector of length n #
        Integral = total/deno # n*MC matrix #

        lgLik = sum(log(deno)) # calculate the log-likelihood #

        post1 = as.matrix(Integral%*%bi)/MC # n*ncz matrix #
        tempB = if(ncz>1) t(apply(bi,1,function(x) x%o%x)) else bi^2 # MC*(ncz^2) matrix #
        post2 = as.matrix(Integral%*%tempB)/MC # n*(ncz^2) matrix #
        BSigma_new = matrix(colMeans(post2),nrow=ncz) # ncz*ncz matrix #

        post3 = rowMeans((Y-Ymu)^2*Integral[ID,])
        Ysigma_new = sum(post3)/N

        tempY = as.vector(rowMeans((Y-Ymu)*Integral[ID,])) # vector of length N #
        post4 = crossprod(tempY,X) # vector of length ncx #

        exp.es = exp(eta.s)
        temp1 = as.vector(rowMeans(exp.es*Integral[Index,])) # vector of length M #
        temp2 = Xtime2*temp1 # M*ncx matrix #
        temp3 = if(ncx>1) t(apply(Xtime2,1,function(x) x%o%x))*temp1 else Xtime2^2*temp1
        # M*(ncx^2) matrix #       
        tempXX = if(ncx>1) t(apply(X,1,function(x) x%o%x)) else X^2
        # N*(ncx^2) matrix # 
        tempXX = if(ncx>1) matrix(colSums(tempXX),nrow=ncx) else sum(tempXX)
        # ncx*ncx symmetric matrix #

        post5 = as.matrix(apply(temp2,2,function(x) tapply(x,Index1,sum))) # n_u*ncx matrix #
        post6 = as.matrix(apply(temp3,2,function(x) tapply(x,Index1,sum))) # n_u*(ncx^2) matrix #

        betaScore = as.vector(post4/Ysigma_new+alpha_old*colSums(d*Xtime)-alpha_old*colSums(lamb_old*post5))
        betaInfo = -tempXX/Ysigma_new-(alpha_old^2)*matrix(colSums(lamb_old*post6),nrow=ncx)
        beta_new = as.vector(beta_old-solve(betaInfo)%*%betaScore)

        eta2n = as.vector(Xtime2%*%beta_new)+Ztime2%*%t(bi) # M*MC matrix #
        eta.sn = as.vector(Wtime%*%gamma_old)[Index]+alpha_old*eta2n # M*MC matrix #
        exp.esn = exp(eta.sn)
        temp4 = as.vector(rowMeans(exp.esn*Integral[Index,])) # vector of length M #
        temp5 = as.vector(rowMeans(eta2n*exp.esn*Integral[Index,])) # vector of length M #
        temp6 = as.vector(rowMeans(eta2n^2*exp.esn*Integral[Index,])) # vector of length M #
        Wtime2 = as.matrix(Wtime[Index,])
        temp7 = Wtime2*temp4 # M*ncw matrix #
        temp8 = Wtime2*temp5 # M*ncw matrix #
        temp9 = if(ncw>1) t(apply(Wtime2,1,function(x) x%o%x))*temp4 else Wtime2^2*temp4
        # M*(ncw^2) matrix # 

        post7 = as.vector(tapply(temp5,Index1,sum)) # vector of length n_u #
        post8 = as.vector(tapply(temp6,Index1,sum)) # vector of length n_u #
        post9 = as.matrix(apply(temp7,2,function(x) tapply(x,Index1,sum))) # n_u*ncw matrix #
        post10 = as.matrix(apply(temp8,2,function(x) tapply(x,Index1,sum))) # n_u*ncw matrix #
        post11 = as.matrix(apply(temp9,2,function(x) tapply(x,Index1,sum))) # n_u*(ncw^2) matrix #

        alphaScore = sum(d*(as.vector(Xtime%*%beta_new)+rowSums(Ztime*post1)))-sum(lamb_old*post7)
        alphaInfo = -sum(lamb_old*post8)
        gammaScore = colSums(d*Wtime)-colSums(lamb_old*post9)
        gammaInfo = -matrix(colSums(lamb_old*post11),nrow=ncw)
        ag_Info = -colSums(lamb_old*post10)
        agScore = c(gammaScore,alphaScore)
        agInfo = matrix(0,nrow=(ncw+1),ncol=(ncw+1))
        agInfo[1:ncw,1:ncw] = gammaInfo
        agInfo[(ncw+1),1:ncw] = ag_Info
        agInfo[1:ncw,(ncw+1)] = ag_Info
        agInfo[(ncw+1),(ncw+1)] = alphaInfo
        ag_old = c(gamma_old,alpha_old)
        ag_new = as.vector(ag_old-solve(agInfo)%*%agScore)
        gamma_new = ag_new[1:ncw]
        alpha_new = ag_new[ncw+1]

        eta.sn2 = as.vector(Wtime%*%gamma_new)[Index]+alpha_new*eta2n # M*MC matrix #
        temp10 = as.vector(rowMeans(exp(eta.sn2)*Integral[Index,])) # vector of length M #
        post12 = as.vector(tapply(temp10,Index1,sum)) # vector of length n_u #
        lamb_new = Index2/post12

        result = list(beta=beta_new,Ysigma=Ysigma_new,BSigma=BSigma_new,gamma=gamma_new,alpha=alpha_new,
                 lamb=lamb_new,lgLik=lgLik)
        
        return(result)
}        



#=============== Initial Value Calculation ===============#

"JMinitECM" <- function(ModelInfo,beta){

         W = ModelInfo$W
         start = ModelInfo$start
         stop = ModelInfo$stop
         event = ModelInfo$event
         Y.est = ModelInfo$Y.est
         data.init = data.frame(start=start,stop=stop,event=event,W=W,Y.est=Y.est)
         fit = coxph(Surv(start,stop,event)~W+Y.est,data=data.init)

         ncw = ncol(W)
         gamma = fit$coefficients[1:ncw]
         alpha = fit$coefficients[ncw+1]
    
         bBLUP = ModelInfo$b
         Xtime2 = ModelInfo$Xtime2
         Wtime = ModelInfo$Wtime
         Index = ModelInfo$Index
         Index1 = ModelInfo$Index1
         Index2 = ModelInfo$Index2
         Ztime2 = ModelInfo$Ztime2
         Wtime2 = as.matrix(Wtime[Index,])
         Ytime2 = as.vector(Xtime2%*%beta)+rowSums(Ztime2*bBLUP[Index,]) # vector of length M #

         temp = as.vector(exp(Wtime2%*%gamma+alpha*Ytime2)) # M*1 vector #
         risk = as.vector(tapply(temp,Index1,sum)) # vector of length n_u #
         lamb = Index2/risk

         result = list(gamma=gamma,alpha=alpha,lamb=lamb)
         return(result)
}




#========== Joint Modeling Main Function Using ECM-algorithm ==========#

# fitMLE      :: Fit from linear mixed effects model
# fitCOX      :: Fit from Cox model
# timeVarY    :: Specs for time-varying variables
# funcY       :: Specs for time-varying variables
# control     :: List of options such as seed (warning: if seed is set then it
#                resets the seed every iteration), tolerances and MC sample size.
# print.every :: How often to print status updates to screen 

"JointModelECM" <- function(fitLME,fitCOX,timeVarY,funcY,control,print.every=1)
{
            if(!inherits(fitLME,"lme")){
               stop("\n'fitLME'must be a fit returned by lme().")
            }
            if(length(fitLME$group)>1){
               stop("\nnested random-effects are not allowed in lme().")
            }
            if(!is.null(fitLME$modelStruct$corStruct)){
               warning("correlation structure in 'fitLME' is ignored.\n")
            }
            if(!is.null(fitLME$modelStruct$varStruct)){
               warning("variance structure in 'fitLME' is ignored.\n")
            }

            if(!inherits(fitCOX,"coxph")){
               stop("\n'fitCOX' must be a fit returned by coxph().")
            }
            if(is.null(fitCOX$x)){
               stop("\nmust specify argument 'x=TRUE' when using coxph().")
            }
            
            ID = as.vector(unclass(fitLME$groups[[1]])) 
            # ids of the subjects, each is repeated by the number of longitudinal measurements #
            ni = as.vector(tapply(ID,ID,length)) 
            # number of longitudinal measurements for each subject #          
            b = data.matrix(ranef(fitLME)) 
            # estimated random effects for each subject #
            dimnames(b) = NULL
            nY = nrow(b)

            if(ncol(fitCOX$y)!=3){
               stop("\nmust fit time-dependent Cox model in coxph().")
            }
            start = as.vector(fitCOX$y[,1])
            stop = as.vector(fitCOX$y[,2])
            event = as.vector(fitCOX$y[,3])
            Time = stop[cumsum(ni)]
            # observed event time, may be right-censored #
            d = event[cumsum(ni)] 
            # censoring index #
            nT = length(Time) 
            # sample size for the survival part #
            if(sum(d)<5)
               warning("\nmore than 5 events are required.")
            if(nY!=nT)
               stop("sample sizes in the longitudinal and event processes differ.\n")

            W = as.matrix(fitCOX$x)
            Wtime = as.matrix(W[cumsum(ni),])
            # design matrix of the covariates in survival part, one row for each subject #


            TermsX = fitLME$terms
            data = fitLME$data[all.vars(TermsX)] # data that's involved in the lme model #
            formYx = formula(fitLME) 
            # the fixed-effects formula used when fitting lme model #
            mfX = model.frame(TermsX,data=data) 
            # when fitting lme, if there's any transformation applied to the response or fixed-effect covariates, #
            # give the transformed data. So it have the same structure as data but the entries may be transformed #
            X = model.matrix(formYx,mfX) 
            # give the fixed-effect covariates used in lme, each covariate has one column, if intercept #
            # term is included, there will be one column of all 1's which corresponds to the intercept #

            formYz = formula(fitLME$modelStruct$reStruct[[1]]) # the random-effect formula used in lme #
            mfZ=model.frame(terms(formYz),data=data)
            # give the covariates that's in the random-effect formula, if only random intercept is used, mfZ is a #
            # data frame with 0 columns and N rows. #
            TermsZ = attr(mfZ,"terms") 
            # give information about the random-effect part, e.g. whether random intercept is included #
            Z = model.matrix(formYz,mfZ) 
            # if random intercept is included, the first column of Z is a column of all 1's, and the remaining #
            # columns are from mfZ. #
            y.long = as.vector(model.response(mfX,"numeric"))
            # give the column in mfX which is considered as response, may be transformed, of length N #

            data.id = data[!duplicated(ID),]
            # pick the first row of each subject in data, nrow=n #

            if(!is.null(timeVarY)){
               if(!all(timeVarY%in%names(data)))
                  stop("\n'timeVarY' does not correspond to columns in the fixed-effect design matrix of 'fitLME'.")
               if(length(timeVarY)!=length(funcY))
                  stop("\n'funcY' should be of the same length with 'timeVarY'.")
               numY = length(funcY)
               data.id[timeVarY] = sapply(1:numY, function(i){funcY[[i]](Time)})
            }

            mfX.id = model.frame(TermsX,data=data.id) 
            # difference with data.id is that some entries may be transformed #
            mfZ.id = model.frame(TermsZ,data=data.id)
            # the columns in data.id which are involved in the random effects part #
            Xtime = model.matrix(formYx,mfX.id) # same structure with X, but with only n rows #
            Ztime = model.matrix(formYz,mfZ.id) # same structure with Z, but with only n rows #
            long = as.vector(c(X%*%fixef(fitLME))+rowSums(Z*b[ID,]))
            # estimated values of longitudinal responses, length N #

            U = sort(unique(Time[d==1])) 
            # ordered uncensored observed event time #
            tempU = lapply(Time,function(t) U[t>=U]) 
            # list with n elements, the ith element indicate the time in U which are less than or equal to Time[i] #
            times = unlist(tempU) # vector of length M #
            nk = sapply(tempU,length) 
            # length of each element in times, vector of length n #
            M=sum(nk)
            Index = rep(1:nY,nk) # repeat 1:n by nk, length M #
            Index0 = match(Time,U)
            # vector of length n, for the uncensored subjects, return the value as nk, otherwise return NA # 
            Index1 = unlist(lapply(nk,seq,from=1)) # vector of length M #
            Index2 = colSums(d*outer(Time,U,"==")) # vector of length n_u #

            data.id2 = data.id[Index,]
            if(!is.null(timeVarY)){
               numY = length(funcY)
               data.id2[timeVarY] = sapply(1:numY, function(i){funcY[[i]](times)})
            }
            mfX2 = model.frame(TermsX, data=data.id2)
            mfZ2 = model.frame(TermsZ, data=data.id2)
            Xtime2 = model.matrix(formYx,mfX2)
            Ztime2 = model.matrix(formYz,mfZ2)

            ModelInfo = list(ID=ID,d=d,b=b,start=start,stop=stop,event=event,Index=Index,Index0=Index0,
                        Index1=Index1,Index2=Index2,U=U,Y=y.long,Y.est=long,X=X,Z=Z,W=W,Xtime=Xtime,
                        Xtime2=Xtime2,Ztime=Ztime,Ztime2=Ztime2,Wtime=Wtime)

            VC = lapply(lapply(fitLME$modelStruct$reStruct, as.matrix),
                        function(x){ 
                          return(x * fitLME$sigma^2)
                          })[[1]]
            
            # the estimated variance-covariance matrix for the random effects  #
            beta = as.vector(fixef(fitLME))
            Ysigma = fitLME$sigma^2 

            surv.init = JMinitECM(ModelInfo,beta)
            gamma = surv.init$gamma
            alpha = surv.init$alpha
            lamb= surv.init$lamb
 
            theta_old = list(beta=beta,Ysigma=Ysigma,BSigma=VC,gamma=gamma,alpha=alpha,lamb=lamb)

            ncz = ncol(Z)
            ncx = ncol(X)
            ncw = ncol(Wtime)
            tol.P = control$tol.P
            tol.L = control$tol.L
            iter = control$max.iter
            MC = control$MC
            seed = control$seed

            Ysigmas = rep(0,iter+1)
            BSigmas = matrix(0,nrow=iter+1,ncol=ncz*ncz)
            betas = matrix(0,nrow=iter+1,ncol=ncx)
            gammas = matrix(0,nrow=iter+1,ncol=ncw)
            alphas = rep(0,iter+1)
            lgLik = rep(0,iter+1)
            err.P = err.L = 1
            step = 1

            while(err.P>tol.P && err.L>tol.L && step<=iter){
              
                  betas[step,] = theta_old$beta
                  Ysigmas[step] = theta_old$Ysigma
                  BSigmas[step,] = c(theta_old$BSigma)
                  gammas[step,] = theta_old$gamma
                  alphas[step] = theta_old$alpha

                  theta_new = JMiterECM1(theta_old,ModelInfo,MC,seed)
                  new.P = c(theta_new$beta,theta_new$Ysigma,theta_new$BSigma,theta_new$gamma,theta_new$alpha)
                  old.P = c(theta_old$beta,theta_old$Ysigma,theta_old$BSigma,theta_old$gamma,theta_old$alpha)
                  err.P = max(abs(new.P-old.P)/(abs(old.P)+tol.P))
                  # add tol.P to avoid zero value of the estimated parameters #

                  lgLik[step+1] = theta_new$lgLik
                  new.L = lgLik[step+1]
                  old.L = lgLik[step]
                  err.L = abs(new.L-old.L)/(abs(old.L)+tol.P)

                  if (step%%print.every == 0){
                        cat(paste("Finished iteration ",step,"... (max. iter=",iter,")\n",sep=""))
                  }

                  step = step+1
                  theta_old = theta_new
            }
            betas[step,] = theta_new$beta
            Ysigmas[step] = theta_new$Ysigma
            BSigmas[step,] = c(theta_new$BSigma)
            gammas[step,] = theta_new$gamma
            alphas[step] = theta_new$alpha

            betas = betas[1:step,]    
            Ysigmas = Ysigmas[1:step]
            BSigmas = BSigmas[1:step,]
            gammas = gammas[1:step,]
            alphas = alphas[1:step] 
            lgLik = lgLik[1:step] 

            thetas = list(betas=betas,Ysigmas=Ysigmas,BSigmas=BSigmas,gammas=gammas,alphas=alphas,lgLik=lgLik)
            result = list(thetas=thetas,theta_hat=theta_new,step=step)
            
            return(result)
}

# Compile for speed? 
library(compiler)
if (do_compile){
  JointModelECM <- cmpfun(JointModelECM)
}


#=============== Simulation Study ===============#

load("Aids.RData")

#========== The case with random intercept and random slope ==========#

# Linear mixed effects model fit for stating values:
fitLME <- lme(sqrt(CD4) ~ obstime+obstime:drug, random = ~ obstime | patient, data = Aids)

# Cox-PH model fit for starting values:
fitCOX <- coxph(Surv(start,stop,event) ~ drug+gender, data = Aids, x = TRUE)

# Set timeVarY and funcY:
timeVarY = "obstime"
funcY = list()
funcY[[1]] = function(x){return(x)}

# Other specs:
control = list(tol.P=0.0001,tol.L=10^(-8),max.iter=200,MC=1000,seed=NULL)

# Fit joint model:
fitJOINT = JointModelECM(fitLME,fitCOX,timeVarY,funcY,control)

# Save results:
save(fitJOINT,file="Result.RData")
