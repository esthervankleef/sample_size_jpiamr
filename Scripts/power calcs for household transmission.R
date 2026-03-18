# Simulations for power calculations for household transmission
require(msm)
require(minqa)
# First need to create a data-frame for n-person household  
# containing fields for 
# state, 
# time...i.e. time observed, 
# intervention   1 or 0
# household id
# For a five person household there are 2^5 =32 possible states (binary representation of state
# indicates who is infected)   


# For msm need to supply qmatrix giving allowed transitions in continuous time
# i.e. gain or loss of one infection - done
#  specify constraints as follows using qconstraints

#  all transitions involving going from 0 to 1 infected are constrained to be equal
#  all transitions involving going from 1 to 2 infected are constrained to be equal
#  all transitions involving going from 2 to 3 infected are constrained to be equal
# etc
#  all transitions involving going from reduction by 1 of number infected are constrained to be equal

# effect of intervention is then to act as multiplier on rate of transmission
# assume it's a constant multiplier

# Can also to simulations using pmatrix.msm to get transition probs from a fitted model
# But suggest we use MatrixExp command and intensity matrix to simulate the data
# noting that  "For a continuous-time homogeneous Markov process with transition 
#intensity matrix Q, the probability of occupying state s at time u + t conditional
#on occupying state r at time u is given by the (r,s) entry of the matrix exp(tQ).

# So functions are needed to  
# i)  construct intensity matrix,Q, for n person household with and without intervention to affect transmission
# should take as input a background rate of acquisition, a per person rate of transmission
# and a multiplicative effect of intervention on transmission.
# output is an intensity matrix


# ii) simulate data at different time points for m households
#   inputs are intensity matrix, initial conditions and time points for simulation, and 
#   number of households  
# output is data frame with household ID, intervention (yes or now), state and time


# iii) 

#  construct qconstraint vector which specifies which baseline transitions intensities are 
# equal for n person household 
# input is n, size of the household
# output is a vector of integers, 1,2,3...

#Noting that 
#"the intensity parameters are assumed to be ordered by reading across the rows of the transition matrix, starting at the first row, ignoring the diagonals.

# Note that for an n person household the state numbers are based on 1+binary representation of states
# where 0 is uninfected and 1 is infected
# so for three person household
# 000  = 0 +1 (state 1, no-one infected)
# 001  = 1 +1 (state 2, person 3 infected)
# 110  = 6 +1 (state 7, person 1 and 2  infected)
# 011  = 3 +1 (state 4, person 2 and 3  infected)
# 



get_num_infected<-function(x, maxhhsize=7){
  # x is state number, hhsize is household size (and max number infected)
  # returns number infected which is given by the number of ones
  # in the binary representation of (x-1)
  y<-as.raw(x-1)
  maxbits<-1+maxhhsize # 
  numberinfected<-0
  for(i in 1:maxbits){
     infected<-as.integer(y & as.raw(2^(i-1)))
     numberinfected<-numberinfected+(infected>0)
  }
  return(numberinfected)
}

infection_added<-function(state1,state2){
  #state1 and state2 are integer representations of two states, 1...n
  # where 1 means none infected, and for >1
  # the state number -1 gives the binary representation of who is infected
  # where 0 is uninfected and 1 is infected
  # function returns true if state 2 can be obtained from state 1
  # by changing one digit from 0 to 1 in the binary representation
  n1<-state1-1
  n2<-state2-1
  # infection is added if n2>n1 and the difference is a power of 2
  if(n2>n1){
     if(log2(n2-n1)%%1 ==0 & bitwAnd(n1,n2)==n1) return(TRUE) else return(FALSE)
  } else {
    return(FALSE)
  } 
}

infection_removed<-function(state1,state2){
  #state1 and state2 are integer representations of two states, 1...n
  # where 1 means none infected, and for >1
  # the state number -1 gives the binary representation of who is infected
  # where 0 is uninfected and 1 is infected
  # function returns true if state 2 can be obtained from state 1
  # by changing one digit from 1 to 0 in the binary representation
  n1<-state1-1
  n2<-state2-1
  # infection is added if n2>n1 and the difference is a power of 2
  if(n2<n1){
    if(log2(n1-n2)%%1 ==0  & bitwAnd(n1,n2)==n2) return(TRUE) else return(FALSE)
  } else {
    return(FALSE)
  } 
}


make_intensity_matrix<-function(hhsize,intervention_effect=1,bg_rate=3,p2p_rate=5,clearance_rate=10){
  # returns intensity matrix,Q, for hhsize person household intervention having multiplicative effect on transmission
  # this matrix to be used for simulation (using matrix exponentiation to determine state transition probabilities)
  # given by intervention_effect 
  #bg_rate is background rate of transmission
  #p2p_rate is rate of transmission from one infected person to another
  # clearance_rate is rate of clearing carriage
  numstates<-2^hhsize
  statenumbers_dec<-1:numstates  # statenumbers in decimal
  statenumbers_raw<-as.raw(statenumbers_dec)  # raw byte
  
  numberinfectedforstate<-get_num_infected(1:numstates, maxhhsize=hhsize) #  number of people infected in given state
  #calculated by counting 1s in binary representation of state number -1
  Q<-matrix(data=0, nrow=numstates,ncol=numstates)
  
  for(i in 1:numstates){
    for(j in 1:numstates){
      if(i!=j){ #deal with diagonal later
        if(infection_removed(i,j)) {
          Q[i,j] <- clearance_rate 
        } else if(infection_added(i,j)) {
          Q[i,j] <- bg_rate +   numberinfectedforstate[i]*p2p_rate*intervention_effect        
        } else{
          # do nothing
        }
      }
    }
    Q[i,i] <- - sum(Q[i,])
  }
  return(Q)
}

samplestate0<-function(Q){
  #choose initial state at random based on long run state of Q
  nstates<-dim(Q)[1]
  longermprobs<-MatrixExp(Q,t=10000)[1,]
  state<-sample(nstates, size=1, prob=longermprobs)
  return(state)
}
  
  
simhhdata<-function(Q,state0, times){
  # inputs Q intensity matrix
  # state0 - state to initalise to 
  # times - sampling times for outputs
  #  returns a vector of sampled dates at times
  nstates<-dim(Q)[1]
  currentstate<-state0
  statesattimes<-NULL
  lasttime<-0
  for(i in 1:length(times)){
      interval<-times[i]-lasttime
      P<-MatrixExp(Q,t=interval)
      transprobs<-P[currentstate,]
      newstate<-sample(nstates, size=1, prob=transprobs)
      statesattimes<-c(statesattimes,newstate)
      lasttime<-times[i]
  }
  return(statesattimes)
}

makefakehhdata<-function(numhh, hhsize,intervention_effect=1,bg_rate=0.5475,p2p_rate=3.869,clearance_rate=0.949,times){
  # numh is number of households, each of size hhsize. Other parameters as above
   n<-length(times)*numhh
    output.df<-data.frame(hhid=rep(NA,n), time=rep(NA,n), treatment=rep(0,n), state=rep(NA,n))
    
    Q<-make_intensity_matrix(hhsize,intervention_effect,bg_rate,p2p_rate,clearance_rate)
    if(intervention_effect!=1) treatment<-TRUE else treatment<-FALSE  
  for(i in 1:numhh){
    nextrow<-(i-1)*length(times)+1
    rows<-nextrow:(nextrow + length(times)-1)
    state0<-samplestate0(Q)
    newstates<-simhhdata(Q, state0,times)
    output.df$hhid[rows]<-i
    output.df$treatment[rows]<- treatment
    output.df$time[rows]<-times
    output.df$state[rows]<-newstates
  }
  return(output.df)
}


makeQconstraintmat<-function(qmatrix){
  # qmatrix is the intensity matrix for the Markov process
#  returns a vector of indicators specifying which baseline transition intensities are equal. For example,
# qconstraint = c(1,2,3,3) constrains the third and fourth intensities to be equal, in a model with four allowed instantaneous transitions.
# do this for covariates set to 0  
  
  # first create these matrixes of same dimensions as qmat
  # 1. matrix m.nonzero that holds true for nonzero elements of qmat, false otherwise
  # 2. matrix m.infadded that holds true if one infection is being added n state transition
  # 3. matrix m.infremoved that holds true if one infection is being removed n state transition
  # 4. matrix m.1to2 that holds true if there is one persons infected and moving to a state with one more infected
  # 5. matrix m.2to3 that holds true if there are two  people infected and moving to a state with one more infected
  # m.constraints matrix of constraints to be added 
  n<-dim(qmatrix)[1]
  hhsize<-log2(n)
  m.nonzero<-qmatrix >0
  m.infadded<-matrix(FALSE, nrow=n, ncol=n)
  m.infremoved<-matrix(FALSE, nrow=n, ncol=n)
  m.0to1<-matrix(FALSE, nrow=n, ncol=n)
  m.1to2<-matrix(FALSE, nrow=n, ncol=n)
  m.2to3<-matrix(FALSE, nrow=n, ncol=n)
  m.3to4<-matrix(FALSE, nrow=n, ncol=n)
  m.4to5<-matrix(FALSE, nrow=n, ncol=n)
  m.5to6<-matrix(FALSE, nrow=n, ncol=n)
  m.6to7<-matrix(FALSE, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      if(infection_removed(i,j)) m.infremoved[i,j]<-TRUE         
      if(infection_added(i,j)){
        switch(get_num_infected(i)+1,
               m.0to1[i,j] <-TRUE,
               m.1to2[i,j] <-TRUE,
               m.2to3[i,j] <-TRUE,
               m.3to4[i,j] <-TRUE,
               m.4to5[i,j] <-TRUE,
               m.5to6[i,j] <-TRUE,
               m.6to7[i,j] <-TRUE
               )
      }
    }
  }
  constraints.matrix<-matrix(NA, nrow=n, ncol=n)
  constraints.matrix[m.0to1] <-1 
  constraints.matrix[m.infremoved] <-2 
  if(hhsize>1)  constraints.matrix[m.1to2] <-3
  if(hhsize>2)  constraints.matrix[m.2to3] <-4
  if(hhsize>3)  constraints.matrix[m.3to4] <-5
  if(hhsize>4)  constraints.matrix[m.4to5] <-6
  if(hhsize>5)  constraints.matrix[m.5to6] <-7
  if(hhsize>6)  constraints.matrix[m.6to7] <-8

  constraints.vector<-t(constraints.matrix)[t(m.nonzero)]
  return(constraints.vector)
}  

# Simulate data 
# Assumptions 12 households per cluster so 264 households in total
# 22 clusters, 11 with the intervention


hhsize<-5
hh<-12*11
# assuming time units of years here

# Baseline parameters - informed by Haverkate et al CMI 2017 https://doi.org/10.1016/j.cmi.2016.08.021
# This paper reports that "Person-to-person transmission occurred at a 
# rate of 0.0053/colonized person/day (0.0025–0.011), background transmission
# at 0.00015/day (95% CI 0.00002–0.00039), and decolonization at 0.0026/day (0.0016–0.0040) for index patients and 0.0090/day (0.0046–0.018) for household members.
# so with timescale of years
# decolonization rate is  365*.0026 = 0.949 /yr 
drate<-0.949
# background transmsission is likely to be substantially higher in sSA than in Netherlands 
# We assume a rate 10 times higher i.e. 0.0015 /day or 0.5475 per year 
brate<-0.5475
# For  person-to-person assume at baseline twice the rate in the Netherlands
# i. 2* 0.0053 = 0.0106 /day (i.e. mean time to transmission of 94 times between 2 people)
# so rate to use is 2*.0053 *365  =  3.869 /yr
p2prate<-3.869
#single simulation
testdata1<-makefakehhdata(hh,hhsize,times=c(.25,0.5,1),intervention_effect=1,bg_rate=brate,p2p_rate=p2prate,clearance_rate=drate)
testdata2<-makefakehhdata(hh,hhsize,times=c(.25,0.5,1),intervention_effect=0.75,bg_rate=brate,p2p_rate=p2prate,clearance_rate=drate)
testdata2$hhid<-testdata2$hhid + max(testdata1$hhid)
testdata<-rbind(testdata1,testdata2)
testdata$num_infected<-get_num_infected(testdata$state)
testdata$logtime<-log(testdata$time)
summary(testdata$num_infected[testdata$treatment==0])
summary(testdata$num_infected[testdata$treatment==1])


require(gee)
fit<-gee(num_infected ~ treatment + logtime , data=testdata, id=hhid, family = poisson, corstr = "AR-M", Mv = 1 )
summary(fit)
temp<-summary(fit)[7]
z<-temp[[1]][2,5]  # this the robust z value
if(z>0) z<- -z
p<-2*(pnorm(z))

get_pvalue_for_intervention<-function(geefit){
  temp<-summary(fit)[7]
  z<-temp[[1]][2,5]  # this the robust z value
  if(z>0) z<- -z
  p<-2*(pnorm(z))
  return(p)
}

#  now repeat above but with many simulations to determine power
numsims<-1000
intervention_effect<- 0.8   # i.e. multiplicative effect on transmission of interventions
hhsize<-5
hh<- 12*11
pvalue_vec<-NULL #vector of pvalues of treatment effect

for(i in 1:numsims){
  print(i)
  testdata1<-makefakehhdata(hh,hhsize,times=c(.25,0.5,1),intervention_effect=1,bg_rate=brate,p2p_rate=p2prate,clearance_rate=drate)
  testdata2<-makefakehhdata(hh,hhsize,times=c(.25,0.5,1),intervention_effect=intervention_effect,bg_rate=brate,p2p_rate=p2prate,clearance_rate=drate)
  testdata2$hhid<-testdata2$hhid + max(testdata1$hhid)
  testdata<-rbind(testdata1,testdata2)
  testdata$num_infected<-get_num_infected(testdata$state)
  testdata$logtime<-log(testdata$time)
  fit<-gee(num_infected ~ treatment + logtime , data=testdata, id=hhid, family = poisson, corstr = "AR-M", Mv = 1 )
  p<-get_pvalue_for_intervention(fit)
  pvalue_vec<-c(pvalue_vec,p)
  
}

# estimated power is proportion of p-values <0.05 under 
power<- sum(pvalue_vec <0.05)/ length(pvalue_vec)
print(c("Estimated power based on ", numsims, " simulations is ", power, "assuming intervention effect ",intervention_effect ), quote=FALSE)





