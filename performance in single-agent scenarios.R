#Scenario 1
true.tox<-c(0.01,0.15,0.45,0.65) # the vector of the true toxicity probabilities for 4 doses
true.eff<-c(0.5,-0.5,-1.5,-3.0) # the vector of the true mean of efficacy responses for 4 doses


prior.tox<-c(0.10,  0.14,  0.18,  0.22) # the vector of prior toxicity probabilities for 4 doses
prior.eff<-c(-1.00, -1.025, -1.05, -1.075) # the vector of prior means of efficacy responses for 4 doses

beta.prior<-rep(3,length(prior.tox)) # the parameters of the Normal-Inverse-Gamma distribution
alpha.prior<-rep(2,length(prior.tox))

set.seed(100)
example<-wdesign(P.tox=true.tox,P.eff=true.eff,
                 target.tox=0.01,target.eff=0.99,n=36,cohort=3,nsims=10000,beta.tox=1,beta.eff=1,beta.prior=beta.prior,alpha.prior=alpha.prior,
                 Pr.tox=prior.tox,Pr.eff=prior.eff,allocation="two.best",delayed.efficacy=F,missing.efficacy=F,
                 safety.constraint=T,rate.tox=0.02,prob.tox.final=0.60,tox.theshold=0.30,prob.tox.theshold=0.95,
                 futulity=T,eff.theshold=0.2,prob.eff.theshold=0.2,rate.eff=0.02,prob.eff.final=0.70,correlation=0.2)
example

##### The proportions of each dose recommendations #####
round(example$recommendations*100,1)

##### The proportions of terminated trials #####
round(example$termination*100,1)

##### The proportions of toxicity and efficacy responses #####
round(example$mean.number.of.tox/example$mean.number.of.patients*100,1)
round(example$mean.number.of.eff/example$mean.number.of.patients,2)
