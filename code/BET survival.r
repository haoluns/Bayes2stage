	
#####################################Code for BET design with survival endpoint#######################################

#####################################Example 1: computing HPD length and posterior probabiilty########################

theta0=1       								     ######### Cutoff for median survival time in stage 1
Pi1=0.8          								 ######### Posterior probability cutoff for stage 1
kHPD=1                                           ######### Value for k_HPD, default set as 1
datat=c(0.1,0.4,0.2,0.6,0.7,0.7,0.1,0.9,0.1,1)   ######### list of event times
is.censored=c(0,1,0,0,0,0,1,0,0,0) 				 ######### list of censoring indicators
		
outputsurvival(Pi1,theta0,kHPD,datat,is.censored)


#####################################Example 2: simulation study for finding average sample size########################

theta0=1         							     ######### Cutoff for median survival time in stage 1
theta1=2         							     ######### Cutoff for median survival time in stage 2
l1=0.5         							         ######### Required HPD interval length for stage 1
l2=0.8           						         ######### Required HPD interval length for stage 2
Pi1=0.7            							     ######### Posterior probability cutoff for stage 1
Pi2=0.8            							     ######### Posterior probability cutoff for stage 2
kHPD=1                                           ######### Value for k_HPD, default set as 1

weibullk=1.5            						 ######### Shape parameter for Weibull distribution
lambda=(theta1+2)^weibullk/log(2)  		         ######### Scale parameter for Weibull distribution, such that the median survival time is theta1+2
censoring=0.8  		                             ######### Assumed rate of censoring
numsim=1000                                      ######### Number of trial replications

simu(censoring,lambda,weibullk,kHPD,theta0,theta1,l1,l2,Pi1,Pi2,numsim)

#####################################Load the following functions into R###############################################

library(R2jags)
library(HDInterval)
library(invgamma)
library(TeachingDemos)



		 
		 
outputsurvival<-function(Pi1,m1,kHPD,datat,is.censored){		
result=weibullposmedian_wrap(datat,is.censored,m1,Pi1)

temp=findparam0(length(datat),Pi1,m1^kHPD/log(2),1-sum(is.censored)/length(datat))

#HPD interval length          
hpdLength = abs((temp[3]*log(2))^(1/kHPD)-(temp[4]*log(2))^(1/kHPD))
#Posterior probability  
posprob = result$prob

return(list(hpd=hpdLength, prob=posprob))
}		 


          
        



hpdresult<-function(meantime,n,conf,censoring){
  a=0.01
  b=0.01
  alpha = a+ censoring * n
  bet = b + n * meantime
  h=hpd(qinvgamma,shape=alpha,rate=bet,conf=conf)
  return(c(h[2]-h[1],h[2],h[1]))
}


poptestvalue<-function(meantime,n,lamda0,censoring){
  
  a=0.01
  b=0.01
  alpha = a+ censoring * n
  bet = b + n * meantime
  
  return(1-pinvgamma(lamda0,shape=alpha,rate=bet))
}







findmeantime<-function(n,conf,lamda0,censoring){
  
  temp<-function(x){poptestvalue(x,n,lamda0,censoring)-conf}
  if (temp(0)*temp(9999)>0){return(999)}
  #if (temp(9999)<0){return(0)}  #add due to weibull error
  
  result = tryCatch({
    y2root=(uniroot(temp,c(0,9999))$root)
  }, warning = function(w) {
    
  }, error = function(e) {
    print(n)
    print(conf)
    print(lamda0)
    print(censoring)
  }, finally = {
    
  })
  
  y2root=(uniroot(temp,c(0,9999))$root)
  return (y2root)
}



findparam0<-function(n,conf,lamda0,censoring){
  
  r=findmeantime(n,conf,lamda0,censoring)
  l=hpdresult(r,n,conf,censoring)
  return(c(r,l))
}




weibullposmedian<-function(datat,is.censored,targetmed,conf=0.9)
{
  
  
  
  N=length(datat)
  
  
  
  params = c( "r","median")
  
  
  status=1-is.censored
  
  t=datat
  is.na(t)<- (status ==0 )
  
  t.cen=datat+status
  
  tinits1<-datat+0.01
  is.na(tinits1)<-status==1
  tinits2<-tinits1+0.01
  
  #inits<-list(list(t=tinits1,r=4,beta=1),list(t=tinits2,r=4,beta=1))
  
  inits<-list(list(t=tinits1,r=1,beta=1))
  data<-list(t=t,t.cen=t.cen,is.censored=is.censored,N=N)
  
  
  
  nc <- 2      #number of MCMC chains to run
  ni <- 10000  #number of samples for each chain     
  nb <- 1000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  
  #testing
  #bugs.out <- jags(data=data, inits=inits, params,n.iter=ni,model.file="mod1.txt", n.chains=1)
  
  
  bugs.out <- jags.model(data=data, inits=inits, file="mod1.txt", n.chains=1)
  codasample=coda.samples(bugs.out,variable.names=params,n.iter=ni)
  #testing 
  
  
  rsample=(codasample[[1]][,2])
  dd=density(rsample)
  restimate=dd$x[which.max(dd$y)]
  #print(bugs.out, digits = 3)
  
  #codasample=coda.samples(bugs.out,variable.names=params,n.iter=ni)
  
  med=(codasample[[1]][,1])
  
  hd=hdi(med,credMass=conf)
  
  
  return(list(sample=med,len=hd[2]-hd[1],k=restimate,prob=sum(med>targetmed)/length(med)))
}


weibullposmedian_wrap<-function(datat,is.censored,target,conf){
  repeatind=1
  countrun=0
  while(repeatind==1 && countrun < 100)
  {
    result <- try(weibullposmedian(datat,is.censored,target,conf));
    if(class(result) != "try-error")
    {  
      repeatind = 0;
    } 
    countrun = countrun + 1 
  }
  if (countrun == 100)return("fail")
  return(result)
}


sink("mod1.txt")        
cat("
    model
{	
    for(j in 1 : N) {                          
    is.censored[j]~dinterval(t[j],t.cen[j])
    t[j]~dweib(r, 1/mu)
    }
    
    beta ~ dgamma(0.01,0.01)
    mu <- 1/beta
    
    median <- pow(log(2) * mu, 1/r)  
    
    r ~ dexp(0.01)
    
}
    ", fill = TRUE)
sink()






#pr(T<C)=0.9
finduniform<-function(lmd,weibullk,target)
{
  #t<-x; c<-y
  #unif=1
  #lmd=1
  #target=0.9
  temp<-function(unif)
  {
    a=	integrate(function(y) { 
      sapply(y, function(y) {
        integrate(function(x) dweibull(x,shape=weibullk,scale=lmd^(1/weibullk))*1/unif, 0, y)$value
      })
    }, 0, unif)$value
    return (a-target)
  }
  
  y2root=(uniroot(temp,c(1,100))$root) 
  return(y2root)
  
}

simu<-function(censoringrate,lambda,weibullk,kHPD,m1,m2,l1,l2,Pi1,Pi2,Numofsim=1000){
  
  
  
  noncensrate=1-censoringrate
  nswitch=99999
  n1=500
  n=1000
  c1=n1*(1-noncensrate)
  c=n*(1-noncensrate)
  n2=n-n1
  c2=c-c1
  
  accuralrate=0.1
  
  runcount=Numofsim
  weibullk=weibullk
  censoringlambda=lambda
  sampleSizeStopped=0
  sampleSizeContinued=0
  
  popsumStage1Lambda0=0
  popsumStage1Lambda1=0
  popsumStage2Lambda0=0
  popsumStage2Lambda1=0
  
  outcomePositive=0
  outcomeNegative=0
  outcomeZero=0
  outcomePositive2=0
  outcomeNegative2=0
  outcomeZero2=0
  
  for (m in 1:Numofsim)
  {
    discard=FALSE
    entrytime = sample(n1)*accuralrate
    
    failuretime=rweibull(n1,shape=weibullk,scale=lambda^(1/weibullk))
    unif=finduniform(lambda,weibullk,1-c1/n1)
    censoringtime=runif(n1,min=0,max=unif)
    
    finalevent=as.integer((censoringtime<failuretime))
    eventtime =  rep(0,n1)
    for (i in 1:length(censoringtime))
    {
      if (finalevent[i]==1)
      {
        eventtime[i]=censoringtime[i]
      }
      else
      {
        eventtime[i]=failuretime[i]
      }
    }
    
    
    realeventtime = entrytime + eventtime
    
    mat=cbind(finalevent,entrytime,eventtime,realeventtime)
    mat=mat[order(mat[,4]),]
    mat=as.data.frame(mat)
    size=length(mat[,1])
    outcome=0
    exittime = 0
    exitsize = 0
    for (i in 1:size)
    {	
      
      currentTime=mat[i,4]
      submat = mat[which(mat[,2]<=currentTime & mat[,4]<=currentTime),] #only focus on number enrolled in the study
      meantime=0
      NumCensor=sum(submat[,1])
      
            
      currentSize=length(submat[,1])
      
      censoring = 1- NumCensor / currentSize
      
      if (sum(1-submat[,1])>=10 && submat[i,1]==0) #only when there is at least 5 death event
      {  
      if (currentSize>0 && censoring >0 )
      {
        datat=submat[,3]
        is.censored=submat[,1]
        
        
        if (sum(1-submat[,1])<nswitch)
        {
          
          
        
        temp=findparam0(currentSize,Pi1,m1^kHPD/log(2),censoring)
        minimumMean = temp[1]
        
        
        hpdLength2 = abs((temp[3]*log(2))^(1/kHPD)-(temp[4]*log(2))^(1/kHPD))
        
        if (hpdLength2<l1)
        {
          result=weibullposmedian_wrap(datat,is.censored,m1,Pi1)
          #result=weibullposmedian(datat,is.censored,m1,Pi1)
          if( result == "fail")
          {
            discard = TRUE
          }else{
          
          #hpdLength = result$len
            hpdLength=hpdLength2  
          poptest =result$prob
          }
        }
        else
        {
          hpdLength=2*l1
          poptest=Pi1
        }
        }
        else if (sum(1-submat[,1])>=nswitch)
        {
          result=weibullposmedian_wrap(datat,is.censored,m1,Pi1)
          
          temp=findparam0(currentSize,Pi1,m1^(result$k)/log(2),censoring)
          minimumMean = temp[1]
          
          
          hpdLength2 = abs((temp[3]*log(2))^(1/result$k)-(temp[4]*log(2))^(1/result$k))
          
            
            #result=weibullposmedian(datat,is.censored,m1,Pi1)
            if( result == "fail")
            {
              discard = TRUE
            }else{
              
              #hpdLength = result$len
              hpdLength=hpdLength2  
              poptest =result$prob
            }
          
          
          
        }
          
      
      if (discard == TRUE)
      {
          runcount = runcount -1 
          outcome=0
          break
      }
      
      if(hpdLength < l1 && poptest > Pi1)
      {
        outcome = 1
        exittime = currentTime
        med = result$sample
        pop1=sum(med>m1)/length(med)
        pop2=sum(med>m2)/length(med)
      exitsize=currentSize
        popsumStage1Lambda0 = popsumStage1Lambda0+pop1
        popsumStage1Lambda1 = popsumStage1Lambda1+ pop2
        break
      }
      else if (hpdLength < l1 && poptest < Pi1)
      {
        outcome = -1
        exittime = currentTime
        exitsize=currentSize
        break
      }else
      {
        
      }
      }
    } } ####stage1 for-loop end
    
    if (outcome == 1 && discard != TRUE )
    {
      outcomePositive = outcomePositive + outcome
      sampleSizeStopped = sampleSizeStopped + exitsize
    } else if (outcome == -1 && discard != TRUE ) {
      outcomeNegative = outcomeNegative + 1
      sampleSizeStopped = sampleSizeStopped + exitsize
    }    else    {
      outcomeZero = outcomeZero + 1
    }
    
    
    if (outcome == 1) #proceed into the next stage
    {
      if (FALSE){
        submatStage1 = mat[which(mat[,2]<=exittime),] # number enrolled in the study for stageone
        
        entrytime2 = sample(n2)*accuralrate + exittime
        
        failuretime2=rweibull(n2,shape=weibullk,scale=lambda^(1/weibullk))
        unif=finduniform(lambda,1-c1/n1)
        censoringtime2=runif(n2,min=0,max=unif)
        
        finalevent2=as.integer((censoringtime2<failuretime2))
        eventtime2 =  rep(0,n2)
        for (i in 1:length(censoringtime2))
        {
          if (finalevent2[i]==1)
          {
            eventtime2[i]=censoringtime2[i]
          }
          else
          {
            eventtime2[i]=failuretime2[i]
          }
        }
        
        
        
        realeventtime2 = entrytime2 + eventtime2
        
        mat2=cbind(finalevent2,entrytime2,eventtime2,realeventtime2)
        colnames(submatStage1)=colnames(mat2)
        mat2=rbind(submatStage1,mat2)
        mat2=mat2[order(mat2[,4]),]
        mat2=as.data.frame(mat2)
        size2=length(mat2[,1])
        outcome2=0
        exittime2 = 0
        
      }#No need for additional enrolling in stage 2 start. we can just use the old one enrolled in stage 1
      
      outcome2=0
      mat2=mat
      size2=length(mat2[,1])
      
      countstart=which(mat2[,4]==exittime)+1
      for (i in countstart:size2)
      {	
        
        currentTime=mat2[i,4]
        submat2 = mat2[which(mat2[,2]<=currentTime & mat2[,4]<=currentTime),] #only focus on number enrolled in the study
        meantime=0
        NumCensor=sum(submat[,1])
        submat=submat2
        
        currentSize=length(submat2[,1])
        censoring = 1- NumCensor / currentSize
        
        if (mat[i,1]==0) #only when there is a death event
        {
        if (currentSize>0  && censoring >0)
        {
          datat=submat[,3]
          is.censored=submat[,1]
          
          if (sum(1-submat[,1])<nswitch)
          {
            
            
            #(findparam0(90,Pi2,m2^kHPD/log(2),censoring)[2]*log(2))^(1/kHPD)
            
            
            temp=findparam0(currentSize,Pi2,m2^kHPD/log(2),censoring)
            minimumMean = temp[1]
            
            
            hpdLength2 =abs((temp[3]*log(2))^(1/kHPD)-(temp[4]*log(2))^(1/kHPD))
          if (hpdLength2<l2)
          {  
            result=weibullposmedian_wrap(datat,is.censored,m2,Pi2)
            
          
            if( result == "fail")
            {
              discard = TRUE
            }else{
            
            #hpdLength = result$len
              hpdLength = hpdLength2
            poptest =result$prob
            }
          } else{
            hpdLength=2*l2
            poptest=Pi2
          }
        
        }
        else if(sum(1-submat[,1])>=nswitch)
        {
          result=weibullposmedian_wrap(datat,is.censored,m2,Pi2)
          
          temp=findparam0(currentSize,Pi2,m2^(result$k)/log(2),censoring)
          minimumMean = temp[1]
          
          
          hpdLength2 = (temp[2]*log(2))^(1/(result$k))
          hpdLength2 = abs((temp[3]*log(2))^(1/result$k)-(temp[4]*log(2))^(1/result$k))
          
          if( result == "fail")
          {
            discard = TRUE
          }else{
            
          hpdLength = hpdLength2
          poptest =result$prob 
          }
        }
        }
        else
        {
          hpdLength=2*l2
          poptest=Pi2
        }
        
        
        if (discard == TRUE)
        {
          runcount=runcount-1
          outcome2 = 0
          exittime2 = currentTime
          exitsize= currentSize
          break
        } 
        
        if(hpdLength < l2 && poptest > Pi2)
        {
          outcome2 = 1
          exittime2 = 
            exitsize= currentSize
          med = result$sample
          pop1=sum(med>m1)/length(med)
        pop2=sum(med>m2)/length(med)
          popsumStage2Lambda0 = popsumStage2Lambda0+pop1
          popsumStage2Lambda1 = popsumStage2Lambda1+pop2
          break
        }
        else if (hpdLength < l2 && poptest < Pi2)
        {
          outcome2 = -1
          exittime2 = currentTime
          exitsize= currentSize
          break
        }else
        {}
        }
      
        
      }###stage 2 for-loop end
      
      if (outcome2 == 1 && discard != TRUE)
      {
        outcomePositive2 = outcomePositive2 + outcome2
        sampleSizeContinued = sampleSizeContinued + exitsize
        
      }
      else if (outcome2 == -1&& discard != TRUE)
      {
        outcomeNegative2 = outcomeNegative2 + 1
        sampleSizeContinued = sampleSizeContinued + exitsize
      }
      else{
        outcomeZero2 = outcomeZero2 + 1
      }
      
    }
  }
  
  return(list(n=sampleSizeContinued/outcomePositive,n1=sampleSizeStopped/runcount))
  
}

