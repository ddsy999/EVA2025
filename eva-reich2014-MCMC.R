
#############################################################:
###    THE MAIN MCMC FUNCTION TO FIT THE HEATWAVE MODEL   ###:
#############################################################:

#############################################################:
#
# INPUTS:
#
#   y      = nyears x ndays matrix of data
#   thresh = the threshold
#   mnB    = prior means (2x5 matrix)
#   sdB    = prior sds (2x5 matrix)
#   initB  = initial values for B
#   temp   = vector of temperatures for heat wave calculations
#   lag    = vector of lags for heat wave calculations
#   iters  = number of MCMC samples to generate
#   burn   = number of discard as burn-in
#   update = sample between graphical displays
#
# OUTPUTS:
#
#   samples = Posterior samples of the model parameters. The first column is the intercept, 
#             the second column is the slope for years.  The five rows correspond to 
#
#             1) probability below the threshold
#             2) GPD scale
#             3) GPD shape
#             4) alpha
#             5) kernel bandwidth
#
#   pwave   = Posterior samples of heatwave probabilities.  pwave[i,j,1] is the ith sample
#             of the probabability of a heatwave of lag[j] days of over temp[j] degrees in
#             year 1.  pwave[i,j,2] is the same for the final year.
#
###################################################################################


prelim2_MCMC<-function(y,thresh,
                       mnB=matrix(0,1,5),sdB=matrix(1,1,5),initB=NULL,
                       temp=100,lag=2,
                       iters=50000,burn=10000,update=100){
  
  # mnB : prior mean of betas 
  # sdB : prior sd of betas 
  # y : n x 5 (x 5) x 4, y[t,(i,)j,r] : (i,)j-th grid, t-th day, r-th run

  # make y as 3-dim array (n x 25 x 4)
  if(length(dim(y)) == 4){
    dim_y <- dim(y)
    y <- array(y, dim = c(dim_y[1], dim_y[2]*dim_y[3], dim_y[4]))
  }
  # Now y[t,k,r] : k-th grid, t-th day, r-th run 
  # y[,1,] : (1,1) grid, y[,2,] : (2,1) grid, ..., y[,6,] : (1,2) grid, ...

    ny <- dim(y)[1] # ny = n is the number of observed days 
    nd <- dim(y)[2] # nd = 25 is the number of grids
    nr <- dim(y)[3] # nr = 4 is the number of runs
    
    # not required 
    # year<-matrix(1:ny,ny,nd,byrow=FALSE)
    # year<-(year-ny/2)/10
    
    if(nd == 5){
      dw2<-as.matrix(dist(1:nd,diag=T,upper=T))^2
    }else if(nd == 25){
      nd_sqrt <- sqrt(nd)
      dw2<-as.matrix(dist(expand.grid(1:nd_sqrt, 1:nd_sqrt),diag=T,upper=T))^2
    }
    

    #INITIAL VALUES:
    nb<-1
    if(is.null(initB)){
      beta<-vector('list', length = 5)
      beta[[1]]<- rep(2, nd)
      beta[[2]]<- rep(1, nd)
      beta[[3]]<-.01
      beta[[4]]<-0
      beta[[5]]<-0
    }


    B<-array(0,c(nd,5)) # making betas; nd (grids) x 5
    # xi, alpha, gamma are the same 
    for(j in 1:5){
      B[,j]<-make.B(beta[[j]],j)
    }
    # prob, sig : varying across location 
    # xi(shape), alpha, gamma : constant 
    
    L <- nd # T_l = 1:nd
    FAC<-array(0,c(nd,L)) # L = nd
    a<-A<-matrix(1,nd) # a : PS distributed B_ell, A : generated random effects A_j
    # for(t in 1:ny){
      FAC<-fac2FAC(make.fac(dw2,B[1,5])) # gaussian kernels  K(j)s
      A<-a2A(FAC,a,B[1,4]) 
    # }

    curll_arr<-array(0,dim = c(ny,nd,nr))
    curlp <- matrix(0, nd)
      curll_arr<-loglike(y,A,B,thresh)
      curll <- apply(curll_arr, 2, sum)
       # current log likelihood
      for(k in 1:nd){
        curlp[k]<-dPS(a[k],0,1,B[1,4]) # current log pdf of random effects 
      }
    

    samples<-array(0,c(iters,length(unlist(beta))))
    dimnames(samples)[[1]]<-paste("Iteration",1:iters)
    
    paste_text <- function(text, grid_mat){
      return(paste0(text, "(", grid_mat[1], ", ", grid_mat[2], ")"))
    }
    
    if(nd == 5){
      dimnames(samples)[[2]]<-c(paste0("Prob", 1:nd),
                                paste0("Scale", 1:nd),
                                "Shape", "Alpha", "BW") 
    }else if(nd == 25){
      grid_mat <- expand.grid(1:nd_sqrt, 1:nd_sqrt)
      dimnames(samples)[[2]]<-c(apply(grid_mat, 1, paste_text, text = "Prob"), 
                                apply(grid_mat, 1, paste_text, text = "Scale"), 
                                "Shape", "Alpha", "BW") 
    }
     
    
    
    
    # nw    <- length(lag)
    # ymin  <- min(year)
    # ymax  <- max(year)
    # pwave <- array(0,c(iters,nw,2))
    # dimnames(pwave)[[1]]<-paste("Iteration",1:iters)
    # dimnames(pwave)[[2]]<-paste("temp",temp,"lag",lag)
    # dimnames(pwave)[[3]]<-c("First year","Last year")


    cuts<-exp(c(-1,0,1,2,5,10))
    MHa<-rep(1,100)
    atta<-acca<-0*MHa
    attb<-accb<-MHb<-matrix(c(.2,.1,.1,.02,.02),nb,5,byrow=T)


    for(i in 1:iters){

      ####################################################
      ##############      Random effects a    ############
      ####################################################

      
      
      olda<-a
      
        parts<-loglikeparts(y,A,B,thresh)
        W<-ifelse(parts$above,1,0) # W : 1 if above, 0 otherwise 
        ccc<--A*apply(parts$expo, 2, sum)+apply(W, 2, sum)*log(A) 
        # if above, ccc = -A*L(y)^(1/alpha) + log(A), otherwise ccc = -A*(log(1/pi))^(1/alpha)
        # summation of log likelihoods 
        alpha<-B[1,4] 


        for(k in 1:L){
         l1<-get.level(a[k],cuts) 
         cana<-exp(rnorm(1,log(a[k]),MHa[l1])) # candidate of B ~ lognormal(mean = log(oldA), scale) 
         # initially scale = 1, changed if accept ratio is too high or low
         l2<-get.level(cana,cuts) 
         WWW<-FAC[,k]^(1/alpha) # weights, K(j)^(1/alpha)
         canA<-A+WWW*(cana-a[k]) # update candidate of A from candidate of B (as A[t,] = WWW * a[t,k])
         cc<--canA*apply(parts$expo, 2, sum)+apply(W, 2, sum)*log(canA)
         canlp<-dPS(cana,0,1,alpha) # density of B (PS), candidate (log) pdf of random effects A  
         R<-sum(cc-ccc)+
            canlp-curlp[k]+ 
            dlognormal(a[k],cana,MHa[l2])-
            dlognormal(cana,a[k],MHa[l1]) # log acceptance ratio, i think 
         if(!is.na(exp(R))){if(runif(1)<exp(R)){ # if accept the new candidate, update... 
             a[k]<-cana 
             A<-canA
             ccc<-cc
             curlp[k]<-canlp
         }}
        }
        
        
    
       curll_arr<-loglike(y,A,B,thresh)
       curll <- apply(curll_arr, 2, sum)
      # update current log likelihood 
      


      ####################################################
      ##############    Model parameters, B   ############
      ####################################################

      
     if(i>25){for(k in 1:5){if(k<3 | j<3){ # nb = 1, so j must be 1? 
       # k < 3 : corresponding parameter is prob(pi) or sig (GPD parameters) 
       # k = 3 : corresponding parameter is xi 
       # k > 3 : corresponding parameter affects dependence structure (alpha, gamma)
       attb[1,k]<-attb[1,k]+1
       canbeta<-beta
       canB<-B
       canA<-A
       canll_arr<-curll_arr
       canlp<-curlp
       canFAC<-FAC
        
       canbeta[[k]]<-beta[[k]]+MHb[1,k]*rnorm(length(beta[[k]])) # candidate of beta? 
       canB[,k]<-make.B(canbeta[[k]],k)
       if(k>3){
         alpha<-canB[1,4]
         gamma<-canB[1,5]
         canFAC<-fac2FAC(make.fac(dw2,gamma))
         canA<-a2A(canFAC,a,alpha) # make new candidate of random effects A from candidate parameters
       }


         canll_arr<-loglike(y,canA,canB,thresh) 
         canll <- apply(canll_arr, 2, sum)# candidate log likelihood 
       }

       if(k==4){for(l in 1:L){
         canlp[l]<-dPS(a[l],0,1,canB[1,4]) # candidate log pdf for alpha (PS)
       }}

       R<-sum(canll-curll)+
          sum(canlp-curlp)+
          sum(dnorm(canbeta[[k]],mnB[1,k],sdB[1,k],log=T))-
          sum(dnorm(beta[[k]],mnB[1,k],sdB[1,k],log=T)) # log accept ratio 

       if(!is.na(exp(R))){if(runif(1)<exp(R)){
          accb[1,k]<-accb[1,k]+1
          beta<-canbeta
          B<-canB
          A<-canA
          curll_arr<-canll_arr
          curll <- canll
          curlp<-canlp
          FAC<-canFAC 
       }}
     }}


     level<-get.level(olda,cuts)
     for(j in 1:length(MHa)){
         acca[j]<-acca[j]+sum(olda[level==j]!=a[level==j])
         atta[j]<-atta[j]+sum(level==j)
         if(i<burn/2 & atta[j]>100){
           if(acca[j]/atta[j]<0.3){MHa[j]<-MHa[j]*0.9}
           if(acca[j]/atta[j]>0.6){MHa[j]<-MHa[j]*1.1}
           acca[j]<-atta[j]<-0
         }
     }


      for(j in 1:nb){for(k in 1:5){if(i<burn/2 & attb[j,k]>50){
        if(accb[j,k]/attb[j,k]<0.3){MHb[j,k]<-MHb[j,k]*0.9}
        if(accb[j,k]/attb[j,k]>0.6){MHb[j,k]<-MHb[j,k]*1.1}
         accb[j,k]<-attb[j,k]<-0
      }}}

      #KEEP TRACK OF STUFF:
      samples[i,]<-unlist(beta)
      # for(j in 1:nw){
      #    pwave[i,j,1]<-heatwave(temp[j],lag[j],B[1,1,],thresh)
      #    pwave[i,j,2]<-heatwave(temp[j],lag[j],B[ny,1,],thresh)
      # }


     #DISPLAY CURRENT VALUE:
     if(i%%update==0){
       print(paste(i, "th update at", Sys.time()))
       par(mfrow=c(5,2),mar=c(2,2,2,2))
       
       if(nd == 5){
         plot(samples[1:i,1], main="Prob 1",type="l")
         lines(samples[1:i,2],main="Prob 2",type="l") 
         lines(samples[1:i,3],main="Prob 3",type="l") 
         lines(samples[1:i,4],main="Prob 4",type="l") 
         lines(samples[1:i,5],main="Prob 5",type="l")
         abline(h=0) 
         
         plot(samples[1:i,5+1], main="Scale 1",type="l")
         lines(samples[1:i,5+2],main="Scale 2",type="l") 
         lines(samples[1:i,5+3],main="Scale 3",type="l") 
         lines(samples[1:i,5+4],main="Scale 4",type="l") 
         lines(samples[1:i,5+5],main="Scale 5",type="l")
         abline(h=0) 
        
         
         plot(samples[1:i,11],main="Shape",type="l") 
         abline(h=0) 
         plot(samples[1:i,12],main="Alpha",type="l")   
         abline(h=0) 
         plot(samples[1:i,13],main="BW",type="l")      
         abline(h=0) 
       }else if(nd == 25){
         plot(samples[1:i,1],main="Prob (1,1)",type="l", ylim = range(samples[1:i, 1:5]), col = "gray70")
         lines(samples[1:i,2],main="Prob (2,1)",type="l", col = "gray70") 
         lines(samples[1:i,3],main="Prob (3,1)",type="l", col = "gray70") 
         lines(samples[1:i,4],main="Prob (4,1)",type="l", col = "gray70") 
         lines(samples[1:i,5],main="Prob (5,1)",type="l", col = "gray70")
         abline(h=0) 
         
         plot(samples[1:i,5],main="Prob (1,5)",type="l", ylim = range(samples[1:i, c(1:5)*5]), col = "gray70")
         lines(samples[1:i,10],main="Prob (2,5)",type="l", col = "gray70") 
         lines(samples[1:i,15],main="Prob (3,5)",type="l", col = "gray70") 
         lines(samples[1:i,20],main="Prob (4,5)",type="l", col = "gray70") 
         lines(samples[1:i,25],main="Prob (5,5)",type="l", col = "gray70")
         abline(h=0) 
         
         plot(samples[1:i,25+1], main="Scale (1,1)",type="l", ylim = range(samples[1:i, 25+1:5]), col = "gray70")
         lines(samples[1:i,25+2],main="Scale (2,1)",type="l", col = "gray70") 
         lines(samples[1:i,25+3],main="Scale (3,1)",type="l", col = "gray70") 
         lines(samples[1:i,25+4],main="Scale (4,1)",type="l", col = "gray70") 
         lines(samples[1:i,25+5],main="Scale (5,1)",type="l", col = "gray70")
         abline(h=0) 
         
         plot(samples[1:i,25+5],  main="Scale (1,5)",type="l", ylim = range(samples[1:i, 25+5*c(1:5)]), col = "gray70")
         lines(samples[1:i,25+10],main="Scale (2,5)",type="l", col = "gray70") 
         lines(samples[1:i,25+15],main="Scale (3,5)",type="l", col = "gray70") 
         lines(samples[1:i,25+20],main="Scale (4,5)",type="l", col = "gray70") 
         lines(samples[1:i,25+25],main="Scale (5,5)",type="l", col = "gray70")
         abline(h=0) 
         
         plot(samples[1:i,51],main="Shape",type="l") 
         abline(h=0) 
         plot(samples[1:i,52],main="Alpha",type="l")   
         abline(h=0) 
         plot(samples[1:i,53],main="BW",type="l")      
         abline(h=0) 
       }
       
       
     }

    }

list(samples=samples)}




#############################################################:
###            OTHER FUNCTION USED IN THE MCMC            ###:
#############################################################:


expit<-function(x){
  1/(1+exp(-x))
}

make.prob<-function(beta){
  # generate pi from beta 
   eta<-beta
   eta[eta>10]<-10
return(expit(eta))}

make.scale<-function(beta){
  # generate sigma from beta 
   eta<-beta
   eta[eta>10]<-10
return(exp(eta))}

make.shape<-function(beta){
  # generate xi from beta 
   eta<-beta
return(eta)}

make.alpha<-function(beta){
  # generate alpha from beta 
   eta<-beta
   eta[eta>10]<-10
return(expit(eta))}

make.gamma<-function(beta){
  # generate gamma from beta 
   eta<-beta
   eta[eta>10]<-10
return(eta)}

make.B<-function(beta,type){
  # generate parameter (of each type) from beta 
  B<-NA
  if(type==1){
    B<-make.prob(beta)
  }
  if(type==2){
    B<-make.scale(beta)
  }
  if(type==3){
    B<-make.shape(beta)
  }
  if(type==4){
    B<-make.alpha(beta)
  }
  if(type==5){
    B<-make.gamma(beta)
  }
return(B)} 

loglikeparts<-function(y,A,B,thresh){  ####TODO

  junk<-is.na(y)
  y<-ifelse(junk,thresh,y)
  
  
  
  if(!is.matrix(B)){
    B <- matrix(B, nrow = 1)
  }
  
   
  if(is.matrix(B)){
    prob<-B[,1]
    sig<-B[,2]
    xi<-B[,3]
    alpha<-B[,4]
  }


  ai<-1/alpha
  above<-y>thresh # TRUE if excess threshold 
  above<-ifelse(junk,FALSE,above) # NA is regarded as censored
  L1<-L3<-expo<-log(1/prob)^ai # (log(1/pi))^(1/alpha)
  
  if(sum(above)>0){
    ttt<-xi*(y-thresh)/sig+1 
    ttt<-ifelse(above,ttt,1) 
    t<-(ttt)^(-1/xi) # s(y) if excess threshold, 1 otherwise 
    l<-1-(1-prob)*t # l(y) or pi_{it}
    # l<-ifelse(above,l,0.5)  
    toohigh<-xi<0 & y>thresh-sig/xi
    t[toohigh]<-l[toohigh]<-ttt[toohigh]<-.1
    L<-log(1/l)
    L2<-L^ai # L^(1/alpha) : exponent 
    L3<-(ai-1)*log(L)+
        log(1-prob)+
        (1+xi)*log(t)-
        log(alpha)-
        log(l)-
        log(sig) # derivative of the exponent 
  }          
  expo<-ifelse(above,L2,L1) # from this l <- ifelse(above, l, 0.5) seems not be necessary 
  expo<-ifelse(junk,0,expo)
  
  log<-ifelse(junk,0,L3)

list(expo=expo,log=log,above=above)}





loglike<-function(y,A,B,thresh){
  
  ny <- dim(y)[1] # ny = n is the number of observed days 
  nd <- dim(y)[2] # nd = 25 is the number of grids
  nr <- dim(y)[3] # nr = 4 is the number of runs

  junk<-is.na(y)
  y<-ifelse(junk,thresh,y)
  
  if(!is.matrix(B)){
    B <- matrix(B, nrow = 1)
  }

  if(is.matrix(B)){
    prob<-B[,1]
    sig<-B[,2]
    xi<-B[,3]
    alpha<-B[,4]
  }

  ai<-1/alpha
  above<-y>thresh
  above<-ifelse(junk,FALSE,above)
  
  
  A_expand <- array(rep(rep(A, each = ny), times = nr), 
                    dim = c(ny,nd,nr))
  B_expand <- array(rep(rep(B, each = ny), times = nr), 
                    dim = c(ny,nd,nr))
  prob_expand <- array(rep(rep(prob, each = ny), times = nr), 
                       dim = c(ny,nd,nr))
  sig_expand <- array(rep(rep(sig, each = ny), times = nr), 
                      dim = c(ny,nd,nr))
  

  lll<- -A_expand*(log(1/prob_expand)^ai)

  if(sum(above)>0){
    ttt<-xi*(y-thresh)/sig+1
    ttt<-ifelse(above,ttt,1)
    t<-(ttt)^(-1/xi)
    l<-1-(1-prob_expand)*t
    # l<-ifelse(above,l,0.5)
    toohigh<-xi<0 & y>thresh-sig/xi
    t[toohigh]<-l[toohigh]<-ttt[toohigh]<-.1
    L<-log(1/l)

    logpabove<--A_expand*(L^ai)+
               log(A_expand)+
               (ai-1)*log(L)+
               log(1-prob_expand)+
               (1+xi)*log(t)-
               log(alpha)-
               log(l)-
               log(sig_expand)
    lll[above]<-logpabove[above]
    lll[toohigh]<--Inf
  }
  lll<-ifelse(junk,0,lll)          
lll}

a2A<-function(FAC,a,alpha){ # a : PS distributed B_ell, A : generated random effects A_j
    #theta is nxnF
    #s is nFxnt
    #alpha in (0,1)
    W<-FAC^(1/alpha)
    if(length(a)==1){xxx<-W*a}
    if(length(a)>1){xxx<-W%*%a}
xxx}  


fac2FAC<-function(x,single=F){
  if(single){x<-x/sum(x)}   
  if(!single){x<-sweep(x,1,rowSums(x),"/")}
x}  

make.fac<-function(dw2,gamma){
   rho2<-exp(gamma)^2
   fac<-exp(-0.5*dw2/rho2)
fac}



rGPD<-function(X,B,thresh){

  if(is.matrix(B)){
    prob<-B[,1]
    sig<-B[,2]
    xi<-B[,3]
  }

  if(!is.matrix(B)){
    prob<-B[1]
    sig<-B[2]
    xi<-B[3]
  }
  
  U<-exp(-1/X)
  U2<-1-(U-prob)/(1-prob)
  Y<-thresh+ifelse(U<prob,0,sig*(U2^(-xi)-1)/xi)
Y}

rGEV<-function(n,mu,sig,xi){
   tau<-runif(n)
   x<--1/log(tau)
   x<-x^(xi)-1
   x<-mu+sig*x/xi
x}




ld<-function(u,A,alpha){
   psi<-pi*u
   c<-(sin(alpha*psi)/sin(psi))^(1/(1-alpha))
   c<-c*sin((1-alpha)*psi)/sin(alpha*psi)
   logd<-log(alpha)-log(1-alpha)-(1/(1-alpha))*log(A)+
         log(c)-c*(1/A^(alpha/(1-alpha)))
exp(logd)}




ld2<-function(u,logs,alpha,shift=0,log=T){

   logs<-logs-shift/alpha
   s<-exp(logs)
   psi<-pi*u
   c<-(sin(alpha*psi)/sin(psi))^(1/(1-alpha))
   c<-c*sin((1-alpha)*psi)/sin(alpha*psi)

   logd<-log(alpha)-log(1-alpha)-(1/(1-alpha))*logs+
         log(c)-c*(1/s^(alpha/(1-alpha)))+
         logs
logd}


rPS<-function(n,alpha,MHA=1,iters=10,initA=NULL){
   
   A<-matrix(0,iters,n)
   accA<-attA<-0
   if(!is.null(initA)){logs<-log(initA)}
   if(is.null(initA)){logs<-rep(0,n)}
   lll<-rep(0,n)
   B<-runif(n)
   lll<-ld2(B,logs,alpha)
   for(i in 1:iters){
     can<-rnorm(n,logs,MHA*(1-alpha))
     ccc<-ld2(B,can,alpha)
     acc<-runif(n)<exp(ccc-lll)
     acc[is.na(acc)]<-F
     logs<-ifelse(acc,can,logs)
     lll<-ifelse(acc,ccc,lll)

     attA<-attA+length(acc)
     accA<-accA+sum(acc)
     if(i<iters/2 & attA>100){
       if(accA/attA<0.4){MHA<-MHA*0.8}
       if(accA/attA>0.5){MHA<-MHA*1.2}
       accA<-attA<-0
     }

     canB<-rtnorm(B)
     ccc<-ld2(canB,logs,alpha)
     R<-ccc-lll+
        dtnorm(B,canB)-
        dtnorm(canB,B)
     acc<-runif(n)<exp(R)
     acc[is.na(acc)]<-F
     B<-ifelse(acc,canB,B)
     lll<-ifelse(acc,ccc,lll)
     A[i,]<-exp(logs)
   }

return(A)}


rtnorm<-function(mn,sd=.2,fudge=0.001){
   upper<-pnorm(1-fudge,mn,sd)
   lower<-pnorm(fudge,mn,sd)
   U<-lower+(upper-lower)*runif(prod(dim(mn)))
return(qnorm(U,mn,sd))}

dtnorm<-function(y,mn,sd=.2,fudge=0.001){
   upper<-pnorm(1-fudge,mn,sd)
   lower<-pnorm(fudge,mn,sd)
   l<-dnorm(y,mn,sd,log=T)-log(upper-lower)
return(l)}


ECkern<-function(h,alpha,gamma,Lmax=50){
    dw2<-rdist(c(0,h),seq(-Lmax,Lmax,1)) 
    W<-fac2FAC(make.fac(dw2,gamma))^(1/alpha)
    for(j in 1:length(h)){
      h[j]<-sum((W[1,]+W[j+1,])^alpha)
    }
h}


npts<-50
Ubeta<-qbeta(seq(0,1,length=npts+1),.5,.5)
MidPoints<-(Ubeta[-1]+Ubeta[-(npts+1)])/2
BinWidth<-Ubeta[-1]-Ubeta[-(npts+1)]

dPS<-function(A,m,s,alpha,npts=100){
   AA<-(A-m)/s
   l<--Inf
   if(AA>0){
      l<-log(sum(BinWidth*ld(MidPoints,AA,alpha)))-log(s)
   }
return(l)}


get.level<-function(A,cuts){
  lev<-A*0+1
  for(j in 1:length(cuts)){
    lev<-ifelse(A>cuts[j],j+1,lev)
  }
return(lev)}

dlognormal<-function(x,mu,sig){
  dnorm(log(x),log(mu),sig,log=T)-log(x)
}

