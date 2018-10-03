library("mvtnorm")


logit<-function(x,alpha.transform=-4.6,beta.transform=-1.5){
  a<-(-1)*beta.transform
  b<-alpha.transform
  y<-exp(-a*x+b)/(1+exp(-a*x+b))
  return(y)
}

loss<-function(p.tox,x.eff,target.tox,target.eff){
  p.eff<-logit(x.eff)
  alpha1<-p.eff*(1-p.tox)
  alpha2<-(1-p.eff)*(1-p.tox)
  alpha3<-1-alpha1-alpha2
  theta1<-target.eff*(1-target.tox)
  theta2<-(1-target.eff)*(1-target.tox)
  theta3<-1-theta1-theta2
  y<-(theta1^2/alpha1)+(theta2^2/alpha2)+(theta3^2/alpha3)-1
  return(y)}


median.beta<-function(x,n){y<-x/n
return(y)}
which.second.max<-function(x){n<-length(x)
y<-which(x==sort(x,partial=n-1)[n-1])
return(y)}
prob.loss<-function(x){h<-length(x)
y1<-1/x
if (any(y1==Inf)){Q<-which(y1==Inf)
y1[Q]<-100}
y2<-mat.or.vec(h,1)
y2[which.max(y1)]<-max(y1)
y2[which.second.max(y1)]<-y1[which.second.max(y1)]
y<-y2/sum(y2)
return(y)}
 
wdesign<-function(P.tox,P.eff,target.tox,target.eff,n,cohort,allocation="one.best",rand.proportion=NULL,missing.efficacy=F,beta.prior,alpha.prior,
                  Pr.tox=NULL,Pr.eff=NULL,initial.prior.value.eff=target.eff,initial.prior.value.tox=target.tox,step.eff=NULL,step.tox=NULL,
                  beta.tox=1,beta.eff=1,kappa=0.5,nsims,safety.constraint=FALSE,toxicity.monotonicity=TRUE,exceptions.list=NULL,rate.tox=NULL,
                  prob.tox.theshold=1,prob.tox.final=NULL,tox.theshold=NULL, coherence=TRUE,coherence.parameter=0,futulity=FALSE,delayed.efficacy=F,
                  efficacy.monotonicity=FALSE,eff.theshold=NULL,prob.eff.theshold=0,rate.eff=NULL,prob.eff.final=NULL,trajectory=F,correlation=0){
  sim<-nsims
  N<-round(n/cohort)                                        
  M<-length(P.tox)
  if (is.null(Pr.tox)){Pr.tox<-c(1:M)*(step.tox/(M-1))+(initial.prior.value.tox-step.tox/(M-1))}    
  if (is.null(Pr.eff)){Pr.eff<-c(1:M)*(step.eff/(M-1))+(initial.prior.value.eff-step.eff/(M-1))}    
  
  counter<-0
  experiment.toxicities<-experiment.efficacies<-experiment<-array(0,dim=c(N+1,M,sim))
  result.exp.tox<-rec.all<-mat.or.vec(sim,M)
  result.tox<-result.eff<-mat.or.vec(sim,1)
  
  exp.tox<-exp.cum.tox<-exp.eff<-exp.cum.eff<-suc.tox<-suc.cum.tox<- suc.eff<- suc.cum.eff<-mat.or.vec(N+1,M)
  eff.outcomes<-mat.or.vec(n+1,M)
  for (z in 1:sim){
    losses<-probability.tox<-probability.eff<-mat.or.vec(N+1,M)
    
    exp.tox[1,]<-exp.cum.tox[1,]<-beta.tox
    suc.tox[1,]<-suc.cum.tox[1,]<-Pr.tox*beta.tox
    exp.eff[1,]<-exp.cum.eff[1,]<-beta.eff
    suc.eff[1,]<-suc.cum.eff[1,]<-Pr.eff*beta.eff
    p.est.tox<-median.beta(suc.cum.tox[1,],exp.cum.tox[1,])
    p.est.eff<-median.beta(suc.cum.eff[1,],exp.cum.eff[1,])
    losses[1,]<-loss(p.est.tox,p.est.eff,target.tox,target.eff)
    
    if(safety.constraint){probability.tox[1,]<-pbeta(tox.theshold, suc.cum.tox[1,]+1,exp.cum.tox[1,]-suc.cum.tox[1,]+1, lower.tail = FALSE)
    for (u in 1:M){if (probability.tox[1,u]>prob.tox.theshold){losses[1,u]<-Inf}}}
    
    # if(futulity){probability.eff[1,]<-pbeta(eff.theshold, suc.cum.eff[1,]+1,exp.cum.eff[1,]-suc.cum.eff[1,]+1, lower.tail = FALSE)
    # check<-prob.eff.theshold   
    # for (u in 1:M){if (probability.eff[1,u]<check){losses[1,u]<-Inf}}}
    
    nextdose<-which.min(losses[1,])
    # cat("nextdose is",nextdose,"\n")
    exp.tox[2,]<-suc.tox[2,]<-0
    exp.tox[2,nextdose]<-cohort
    exp.cum.tox[2,]<-exp.cum.tox[1,]
    exp.cum.tox[2,nextdose]<-exp.cum.tox[1,nextdose]+exp.tox[2,nextdose]
    y<-pnorm(rmvnorm(cohort,mean=c(0,0),sigma=rbind(c(1,correlation),c(correlation,1))))
    toxicities<-(y[,1]<P.tox[nextdose])
    # toxicities<-rbinom(cohort,1,P.tox[nextdose])
    suc.tox[2,nextdose]<-sum(toxicities)
    suc.cum.tox[2,]<-suc.cum.tox[1,] 
    suc.cum.tox[2,nextdose]<-suc.cum.tox[1,nextdose]+suc.tox[2,nextdose]
    
    if(missing.efficacy){rest<-cohort-sum(toxicities)}else{rest<-cohort}
    
    if(!delayed.efficacy){
    exp.eff[2,]<-suc.eff[2,]<-0
    exp.eff[2,nextdose]<-rest
    exp.cum.eff[2,]<-exp.cum.eff[1,]
    exp.cum.eff[2,nextdose]<-exp.cum.eff[1,nextdose]+exp.eff[2,nextdose]
    if(rest>0){
    ind1<-2
      ind2<-2+rest-1
    # eff.outcomes[ind1:ind2,nextdose]<-eff.response<-rnorm(rest,mean=P.eff[nextdose],sd=1)
      if(!missing.efficacy){
eff.outcomes[ind1:ind2,nextdose]<-eff.response<-qnorm(y[,2], mean = P.eff[nextdose], sd = 1)  #no missing
}else{
  eff.outcomes[ind1:ind2,nextdose]<-eff.response<-qnorm(y[which(y[,1]>P.tox[nextdose]),2], mean = P.eff[nextdose], sd = 1)  #missing
      }
    }else{
      ind1<-2
      ind2<-2+cohort-1
      eff.outcomes[ind1:ind2,nextdose]<-eff.response<-rep(0,cohort)
}
      suc.eff[2,nextdose]<-sum(eff.response)
    suc.cum.eff[2,]<-suc.cum.eff[1,] 
    suc.cum.eff[2,nextdose]<-suc.cum.eff[1,nextdose]+suc.eff[2,nextdose]
    }else{
      exp.cum.eff[2,]<-exp.cum.eff[1,]
      suc.cum.eff[2,]<-suc.cum.eff[1,]
    }
    
    previous.dose<-nextdose
    y.previous<-y
    # exp.cum.eff[2,]<-exp.cum.eff[1,]
    # suc.cum.eff[2,]<-suc.cum.eff[1,] 
    j<-2
    
    while (j<N+1){
      prevdose<-nextdose
      p.est.tox<-median.beta(suc.cum.tox[j,],exp.cum.tox[j,])
      p.est.eff<-median.beta(suc.cum.eff[j,],exp.cum.eff[j,])
      losses[j,]<-losses.store<-loss(p.est.tox,p.est.eff,target.tox,target.eff)
      # cat(losses[j,],"\n")
      if(safety.constraint){for (v in 1:M){probability.tox[j,v]<-pbeta(tox.theshold, suc.cum.tox[j,v]+1,exp.cum.tox[j,v]-suc.cum.tox[j,v]+1, lower.tail = FALSE)
      check<-max(prob.tox.theshold-rate.tox*exp.cum.tox[j,v],prob.tox.final)
      if (probability.tox[j,v]>check){losses[j,v]<-Inf}}
        if(toxicity.monotonicity){if(is.null(exceptions.list)){for (v in 1:M){if(losses[j,v]==Inf){if(v!=M){v2<-v+1
        losses[j,v2:M]<-Inf}}}}else{
          for (v in 1:M){
            nochange<-NULL
            for (s in 1:length(exceptions.list)){
              exceptions<-exceptions.list[[s]]
              if(v %in% exceptions){
                nochange<-c(nochange,exceptions)
              }
            }
            exceptions<-unique(nochange)
            if(losses[j,v]==Inf & v!=M){
              change<-seq(v+1,M,1)
              change<-change[!(change %in% exceptions)]
              losses[j,change]<-Inf
            }
          }
        }
        }
      }
      
      if(futulity){
        
        for (v in 1:M){
        vect<-eff.outcomes[,v]  
        vect<-vect[vect!=0]

       if(length(vect)==0){
         mean.vect<-0
         sd.vect<-0
       }else{
         mean.vect<-mean(vect)
         sd.vect<-sum((vect-mean.vect)^2)
       }
        dist.mean<-p.est.eff[v]
        beta0<-beta.prior[v]+sd.vect/2+(exp.cum.eff[j,v]-beta.eff)*beta.eff/exp.cum.eff[j,v]/2*(mean.vect-Pr.eff[v])^2
        lambda0<-exp.cum.eff[j,v]
        alpha0<-alpha.prior[v]+(exp.cum.eff[j,v]-beta.eff)/2
        vari<-beta0/(lambda0*(alpha0-1))
        probability.eff[j,v]<-pnorm(eff.theshold, mean = dist.mean, sd = sqrt(vari), lower.tail = T)
        # cat("dist.mean=",dist.mean,"\n")
        # cat("vari=",vari,"\n")
        # cat("prob=",probability.eff[j,v],"\n")
        
        check<-min(prob.eff.theshold+rate.eff*exp.cum.eff[j,v],prob.eff.final)
        if(probability.eff[j,v]<check){
          losses[j,v]<-Inf
        }
        }
        # cat("probabilities=",probability.eff[j,],"\n")
      }
      
      
      
      if (all(losses[j,]==Inf)){break}  
      
      loss.cand<-losses[j,]
      
      if(toxicity.monotonicity & coherence){
        if(is.null(exceptions.list)){if(sum(toxicities)>coherence.parameter){if(prevdose!=M){prevdose2<-prevdose+1
      loss.cand[prevdose2:M]<-Inf}}else{if(prevdose!=1){prevdose2<-prevdose-1
      loss.cand[1:prevdose2]<-Inf}}}else{
        nochange<-NULL
        for (s in 1:length(exceptions.list)){
          exceptions<-exceptions.list[[s]]
          if(prevdose %in% exceptions){
            nochange<-c(nochange,exceptions)
          }
        }
        exceptions<-unique(nochange)
        if(sum(toxicities)>coherence.parameter){if(prevdose!=M){prevdose2<-prevdose+1
        loss.cand[prevdose2:M]<-Inf
        if(prevdose %in% exceptions){loss.cand[exceptions]<-losses[j,exceptions]}}}else{
          if(prevdose!=1){prevdose2<-prevdose-1
          loss.cand[1:prevdose2]<-Inf
          if(prevdose %in% exceptions){loss.cand[exceptions]<-losses[j,exceptions]}}}
        
      }
      }
      
      #       if (all(loss.cand==Inf)){
      #         losses[j,]<-Inf
      #         break} 
      
      if(allocation=="one.best"){if (all(loss.cand==Inf)){nextdose<-prevdose}else{nextdose<-which.min(loss.cand)}}
      if(allocation=="two.best"){if (all(loss.cand==Inf)){nextdose<-prevdose}else{nextdose<-sample(1:M, 1, prob = prob.loss(loss.cand),replace = TRUE)}}
      if(allocation=="hybrid.1"){N2<-rand.proportion*N
      if(j<N2){nextdose<-which.min(loss.cand)}else{nextdose<-sample(1:M, 1, prob = prob.loss(loss.cand),replace = TRUE)}}
      
      exp.tox[j+1,]<-suc.tox[j+1,]<-0
      exp.tox[j+1,nextdose]<-cohort
      exp.cum.tox[j+1,]<-exp.cum.tox[j,]
      exp.cum.tox[j+1,nextdose]<-exp.cum.tox[j,nextdose]+exp.tox[j+1,nextdose]
      y<-pnorm(rmvnorm(cohort,mean=c(0,0),sigma=rbind(c(1,correlation),c(correlation,1))))
      toxicities.new<-(y[,1]<P.tox[nextdose])
      # toxicities.new<-rbinom(cohort,1,P.tox[nextdose])
      suc.tox[j+1,nextdose]<-sum(toxicities.new)
      suc.cum.tox[j+1,]<-suc.cum.tox[j,] 
      suc.cum.tox[j+1,nextdose]<-suc.cum.tox[j,nextdose]+suc.tox[j+1,nextdose]
      
      #no delayed efficacy
      
      if(delayed.efficacy){
        if(missing.efficacy){rest<-cohort-sum(toxicities)}else{rest<-cohort}
        
        exp.eff[j,]<- suc.eff[j,]<-0
        exp.eff[j,previous.dose]<-rest
        exp.cum.eff[j,]<-exp.cum.eff[j-1,]
        exp.cum.eff[j,previous.dose]<-exp.cum.eff[j-1,previous.dose]+exp.eff[j,previous.dose]
        if(rest>0){
          ind1<-cohort*(j-1)-1
          ind2<-cohort*(j-1)-1+rest-1
          # eff.outcomes[ind1:ind2,nextdose]<-eff.response<-rnorm(rest,mean=P.eff[nextdose],sd=1)
          if(!missing.efficacy){
            eff.outcomes[ind1:ind2,previous.dose]<-eff.response<-qnorm(y.previous[,2], mean = P.eff[previous.dose], sd = 1) # no missing
          }else{
            eff.outcomes[ind1:ind2,previous.dose]<-eff.response<-qnorm(y.previous[which(y.previous[,1]>P.tox[previous.dose]),2], mean = P.eff[previous.dose], sd = 1) # no missing
          }
          suc.eff[j,previous.dose]<-sum(eff.response)
        }
        else{
          ind1<-cohort*(j-1)-1
          ind2<-cohort*(j-1)-1+cohort-1
          eff.outcomes[ind1:ind2,previous.dose]<-eff.response<-rep(0,cohort)
          suc.eff[j,previous.dose]<-sum(eff.response)
          }
        suc.cum.eff[j,]<-suc.cum.eff[j-1,] 
        suc.cum.eff[j,previous.dose]<-suc.cum.eff[j-1,previous.dose]+suc.eff[j,previous.dose]    
        
        exp.cum.eff[j+1,]<-exp.cum.eff[j,]
        suc.cum.eff[j+1,]<-suc.cum.eff[j,] 
      }else{
        if(missing.efficacy){rest<-cohort-sum(toxicities.new)}else{rest<-cohort}
        exp.eff[j+1,]<- suc.eff[j+1,]<-0
        exp.eff[j+1,nextdose]<-rest
        exp.cum.eff[j+1,]<-exp.cum.eff[j,]
        exp.cum.eff[j+1,nextdose]<-exp.cum.eff[j,nextdose]+exp.eff[j+1,nextdose]
        if(rest>0){
        ind1<-cohort*j-1
        ind2<-cohort*j-1+rest-1
        # eff.outcomes[ind1:ind2,nextdose]<-eff.response<-rnorm(rest,mean=P.eff[nextdose],sd=1)
        if(!missing.efficacy){
        eff.outcomes[ind1:ind2,nextdose]<-eff.response<-qnorm(y[,2], mean = P.eff[nextdose], sd = 1) # no missing
        }else{
eff.outcomes[ind1:ind2,nextdose]<-eff.response<-qnorm(y[which(y[,1]>P.tox[nextdose]),2], mean = P.eff[nextdose], sd = 1) # no missing
      }
        }else{
          ind1<-cohort*j-1
          ind2<-cohort*j-1+cohort-1
          eff.outcomes[ind1:ind2,nextdose]<-eff.response<-rep(0,cohort)
        }
        
        suc.eff[j+1,nextdose]<-sum(eff.response)
        suc.cum.eff[j+1,]<-suc.cum.eff[j,] 
        suc.cum.eff[j+1,nextdose]<-suc.cum.eff[j,nextdose]+suc.eff[j+1,nextdose]    
        
      }
      
      # 
      # 
      # exp.eff[j+1,]<- suc.eff[j+1,]<-0
      # exp.eff[j+1,nextdose]<-rest
      # exp.cum.eff[j+1,]<-exp.cum.eff[j,]
      # exp.cum.eff[j+1,nextdose]<-exp.cum.eff[j,nextdose]+exp.eff[j+1,nextdose]
      # ind1<-cohort*j-1
      # ind2<-cohort*j-1+rest-1
      # eff.outcomes[ind1:ind2,nextdose]<-eff.response<-rnorm(rest,mean=P.eff[nextdose],sd=1)
      # suc.eff[j+1,nextdose]<-sum(eff.response)
      # # suc.eff[j,previous.dose]<-sum(rbinom(rest,1,P.eff[previous.dose]))
      # 
      # suc.cum.eff[j+1,]<-suc.cum.eff[j,] 
      # suc.cum.eff[j+1,nextdose]<-suc.cum.eff[j,nextdose]+suc.eff[j+1,nextdose]    
      # 
      # exp.cum.eff[j+1,]<-exp.cum.eff[j,]
      # suc.cum.eff[j+1,]<-suc.cum.eff[j,] 
      toxicities<-toxicities.new
      previous.dose<-nextdose
      y.previous<-y
      j<-j+1}
    
    if (all(losses[j,]==Inf)){counter<-counter+1
    rec.all[z,]<-0
    result.exp.tox[z,]<-(exp.cum.tox[j,]-exp.cum.tox[1,])
    result.tox[z]<-sum(suc.cum.tox[j,]-suc.cum.tox[1,])
    result.eff[z]<-sum(suc.cum.eff[j,]-suc.cum.eff[1,])
    } else{result.exp.tox[z,]<-(exp.cum.tox[N+1,]-exp.cum.tox[1,])
    result.tox[z]<-sum(suc.cum.tox[N+1,]-suc.cum.tox[1,])
    result.eff[z]<-sum(suc.cum.eff[N+1,]-suc.cum.eff[1,])
    j<-N+1
    
    if(delayed.efficacy){
      if(missing.efficacy){rest<-cohort-sum(toxicities)}else{rest<-cohort}
      
      exp.eff[j,]<- suc.eff[j,]<-0
      exp.eff[j,previous.dose]<-rest
      exp.cum.eff[j,]<-exp.cum.eff[j-1,]
      exp.cum.eff[j,previous.dose]<-exp.cum.eff[j-1,previous.dose]+exp.eff[j,previous.dose]
      if(rest>0){
        ind1<-cohort*(j-1)-1
        ind2<-cohort*(j-1)-1+rest-1
        # eff.outcomes[ind1:ind2,nextdose]<-eff.response<-rnorm(rest,mean=P.eff[nextdose],sd=1)
        if(!missing.efficacy){
          eff.outcomes[ind1:ind2,previous.dose]<-eff.response<-qnorm(y.previous[,2], mean = P.eff[previous.dose], sd = 1) # no missing
        }else{
          eff.outcomes[ind1:ind2,previous.dose]<-eff.response<-qnorm(y.previous[which(y.previous[,1]>P.tox[previous.dose]),2], mean = P.eff[previous.dose], sd = 1) # no missing
        }
        suc.eff[j,previous.dose]<-sum(eff.response)
      }
      else{
        ind1<-cohort*(j-1)-1
        ind2<-cohort*(j-1)-1+cohort-1
        eff.outcomes[ind1:ind2,previous.dose]<-eff.response<-rep(0,cohort)
        suc.eff[j,previous.dose]<-sum(eff.response)
        }
      suc.cum.eff[j,]<-suc.cum.eff[j-1,] 
      suc.cum.eff[j,previous.dose]<-suc.cum.eff[j-1,previous.dose]+suc.eff[j,previous.dose]   
    }
    
    p.est.tox<-median.beta(suc.cum.tox[j,],exp.cum.tox[j,])
    p.est.eff<-median.beta(suc.cum.eff[j,],exp.cum.eff[j,])
    losses[j,]<-loss(p.est.tox,p.est.eff,target.tox,target.eff)
    
    if(safety.constraint){for (v in 1:M){probability.tox[j,v]<-pbeta(tox.theshold, suc.cum.tox[j,v]+1,exp.cum.tox[j,v]-suc.cum.tox[j,v]+1, lower.tail = FALSE)
    # check<-max(prob.tox.theshold-rate.tox*exp.cum.tox[j,v],prob.tox.final)
    check<-prob.tox.final

    if (probability.tox[j,v]>check){losses[j,v]<-Inf}}
      if(toxicity.monotonicity){if(is.null(exceptions.list)){for (v in 1:M){if(losses[j,v]==Inf){if(v!=M){v2<-v+1
      losses[j,v2:M]<-Inf}}}}else{
        for (v in 1:M){
          nochange<-NULL
          for (s in 1:length(exceptions.list)){
            exceptions<-exceptions.list[[s]]
            if(v %in% exceptions){
              nochange<-c(nochange,exceptions)
            }
          }
          exceptions<-unique(nochange)
          if(losses[j,v]==Inf & v!=M){
            change<-seq(v+1,M,1)
            change<-change[!(change %in% exceptions)]
            losses[j,change]<-Inf
          }
        }
      }
      }
    }
    
    # if(futulity){for (v in 1:M){probability.eff[j,v]<-pbeta(eff.theshold, suc.cum.eff[j,v]+1,exp.cum.eff[j,v]-suc.cum.eff[j,v]+1, lower.tail = FALSE)
    # check<-min(prob.eff.theshold+rate.eff*exp.cum.eff[j,v],prob.eff.final)
    # if (probability.eff[j,v]<check){losses[j,v]<-99}}
    #   if(efficacy.monotonicity){for (v in 1:M){if(losses[j,v]==99){if(v==1){losses[j,v]<-Inf}else{losses[j,1:v]<-Inf}}}}else{
    #     for (v in 1:M){if(losses[j,v]==99){losses[j,v]<-Inf}}}}
    
    if(futulity){
      
      for (v in 1:M){
        vect<-eff.outcomes[,v]  
        vect<-vect[vect!=0]
        
        if(length(vect)==0){
          mean.vect<-0
          sd.vect<-0
        }else{
          mean.vect<-mean(vect)
          sd.vect<-sum((vect-mean.vect)^2)
        }
        dist.mean<-p.est.eff[v]
        beta0<-beta.prior[v]+sd.vect/2+(exp.cum.eff[j,v]-beta.eff)*beta.eff/exp.cum.eff[j,v]/2*(mean.vect-Pr.eff[v])^2
        lambda0<-exp.cum.eff[j,v]
        alpha0<-alpha.prior[v]+(exp.cum.eff[j,v]-beta.eff)/2
        vari<-beta0/(lambda0*(alpha0-1))
        probability.eff[j,v]<-pnorm(eff.theshold, mean = dist.mean, sd = sqrt(vari), lower.tail = T)
        check<-prob.eff.final
        if(probability.eff[j,v]<check){
          losses[j,v]<-Inf
        }
      }
      # cat("probabilities=",probability.eff[j,],"END","\n")
    }
    
    if (all(losses[j,]==Inf)){counter<-counter+1
    rec.all[z,]<-0
    result.exp.tox[z,]<-(exp.cum.tox[j,]-exp.cum.tox[1,])
    result.tox[z]<-sum(suc.cum.tox[j,]-suc.cum.tox[1,])
    result.eff[z]<-sum(suc.cum.eff[j,]-suc.cum.eff[1,])}else{nextdose<-which.min(losses[j,])
    rec.all[z,nextdose]<-1
    experiment[,,z]<-exp.tox
    experiment.toxicities[,,z]<-suc.tox
    experiment.efficacies[,,z]<-suc.eff}}}
  y<-colSums(rec.all)/sim
  
  if(trajectory){output<-list(all.outcomes=experiment,all.toxicities=experiment.toxicities,all.efficacies=experiment.efficacies,number.of.simulation=sim,True.Toxicity=P.tox,True.Efficacy=P.eff,experimentation=colSums(result.exp.tox)/(n*sim),recommendations=y,termination=counter/sim,mean.number.of.patients=mean(rowSums(result.exp.tox)),mean.number.of.tox=mean(result.tox),mean.number.of.eff=mean(result.eff))}else{
    output<-list(number.of.simulation=sim,True.Toxicity=P.tox,True.Efficacy=P.eff,experimentation=colSums(result.exp.tox)/(n*sim),recommendations=y,termination=counter/sim,mean.number.of.patients=mean(rowSums(result.exp.tox)),mean.number.of.tox=mean(result.tox),mean.number.of.eff=mean(result.eff))}
  return(output)}


##### Example #####
# set.seed(100)
# example<-wdesign(P.tox=c(0.05,0.175,0.20,0.40),P.eff=c(0.10,0.40,0.80,0.95),
#                  target.tox=0.01,target.eff=0.99,n=20,cohort=2,nsims=10000,
#                  Pr.tox=c(0.10,0.20,0.30,0.40),Pr.eff=c(0.65,0.70,0.75,0.80))
# example
