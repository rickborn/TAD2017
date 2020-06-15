#        R Code for Colada 88 - Hot hand artifact
#
# Written by Uri Simonsohn (urisohn@gmail.com)
# date: 2020 05 08
#
#
# This script simulates a sequence of shot attempts, where the probability of conversion
# depends on whether the previous 0,1,2,3 shots were made, and then the correlation between 
# consecutive shots is computed.
# It is used as calibration for the  high correlation between previous shot and prediction
# of next shot in Gilovich Vallone and Tversky Study 4.
#
# Two alterantive definitios of being hot: a streak of 1,2 or 3, whereas the streak is killed and is cold with one miss 
# or the simple average conversion of the last 3 shots. Both lead to the same result
#
###############################################################################################


#Clear everything
  rm(list = ls())

#Function 1 - Compute how long a streak over the last 3 attempts
        get.streak=function(h0)  #h0 is a sequence of 1/0s so far, e.g., h0=c(0,1) means missed 1st, made 2nd shot
        {
        #How many shots so far: h0 is the sequence so far
          k0=length(h0)          
        #If the sequence is in one of the initial 3 shots, we add missed shots before, so players start cold. 
		#Add 0s so at least 3 in the sequence
          h0=c(rep(0,max(3-k0,0)),h0)  # c(1)-->c(0,0,1)
        
        #Adjust k, the number of shots total so far, as we now know it is at least 3
          k=max(k0,3)
        
        #Count consecutive successes prior to this shot
          if (h0[k]==0) return(0)                            #if missed last, streak is 0
          if (h0[k]==1 &  h0[k-1]==0) return(1)              #if made exactly 1 and missed previous, streak is 1
          if (h0[k]==1 &  h0[k-1]==1 & h0[k-2]==0) return(2) #if made previous 2, ...2
          if (h0[k]==1 &  h0[k-1]==1 & h0[k-2]==1) return(3) #if made previous 3 ...3
      }
      
      #Examples of how the function works
          get.streak(c(1,0,1)) #Made, missed, made--> streak of 1 prior to this shot 
          get.streak(c(0,1,1)) #2 in a row
          get.streak(c(1,1,0)) #0 in a row before next throw
          
          
#Function 2 - get.moving.sum - This is an alterantive definition of hand-hotness 
#Instead of setting the streak to 0 if one is missed, the hand is assumed to cool off incrementaly
#so here we compute the moving average of successes over the previous 3 shots      
		get.moving.sum=function(h0)  #h0 is a sequence of 1/0s so far, e.g., h0=c(0,1) means missed 1st, made 2nd shot
        {
        #How many shots so far
          k0=length(h0)
        #Add 0s so at least 3 in the sequence
          h0=c(rep(0,max(3-k0,0)),h0)  # c(1)-->c(0,0,1)
        
        #Adjust k as we now know it is at least 3
          k=max(k0,3)
        
        #Sum last 3 shots
          moving.sum=(h0[k]+h0[k-1]+h0[k-2])
          return(moving.sum)
      }
          
      #Examples
          get.moving.sum(c(1,0,1)) #2 out of 3
          get.moving.sum(c(0,1,1)) #2 out of 3
          get.moving.sum(c(1,1,0)) #2 out of 3
          get.moving.sum(c(1,1,1,1)) #3 out of 3
          get.moving.sum(c(1,1,1,0)) #2 out of 3

                
          
#Function 3 - get the shooting probability for the next shot based on the previous 3 (here we generate the hand-hotness)
		get.shooting.prob=function(h0, p0,p1,p2,p3,type="streak") {
				#SYNTAX
				#  p0,p1,p2,p3 are the probabilities of making the next shot going from the coldest state, 0, to hottest 3
				#  type=c(streak, moving.sum) we define how we consider the hot hand
				
		#STRICT STREAK 
			if (type=="streak")
			{
			streak=get.streak(h0)
			if (streak==0) return(p0)
			if (streak==1) return(p1)
			if (streak==2) return(p2)
			if (streak==3) return(p3)
		  }
		#MOVING AVERAGE
		if (type=='moving.sum')
			{
			moving.sum=get.moving.sum(h0)
			if (moving.sum==0) return(p0)
			if (moving.sum==1) return(p1)
			if (moving.sum==2) return(p2)
			if (moving.sum==3) return(p3)
			}#End if type='moving_average'
	}#End function get shooting ptob
  
  
    #Example
      get.shooting.prob(c(1,0,1),.1,.2,.3,.4,type='streak')
      get.shooting.prob(c(1,0,1),.1,.2,.3,.4,type='moving.sum')
      get.shooting.prob(c(1,0,0),.1,.2,.3,.4,type='streak')
      get.shooting.prob(c(1,1,1),.1,.2,.3,.4,type='streak')
      get.shooting.prob(c(0,1,1),.1,.2,.3,.4,type='streak')

            
#Function 4 - generate sequence of hits/misses, defining how many attempts total, and the probabilities from cold to hot, p0-p4, and how the hotness is defined           
  gen.shots=function(N,p0,p1,p2,p3,type)  {
    h=c()
    for (k in 1:N) h[k]=rbinom(1,size=1,prob=get.shooting.prob(h,p0=p0,p1=p1,p2=p2,p3=p3))
    return(h)
    }

  
#Function 5 - Computes the correlation between a vector and the previous element of the same vector
  cor.real=function(shots)
    {
    shots.previous=shots[2:length(shots)]
    cor.test(shots.previous,shots[1:length(shots)-1])
  }
  
#Function 6 simulate the shots, compute the correlation and mean hit rate, and report results 
	sim=function(...)
	{
    shots=gen.shots(...)
    r=cor.real(shots)$estimate
    p.all=mean(shots)
    return(data.frame(r=as.numeric(r),mean=mean(shots)))
    }

#######################################################
 
 #Calibrations
	#Result in the post 
    #After 0 - 'top' 129   .376 Tarean Prince
    #After 1 - top    75   .450 Bobby Portis
    #After 2 - top    25   .508 Karl-Anthony towns
    #After 3 - Best        .742 Mitchell Robinson
	
   #Method to define hot-hand: sum of hits (hotness is base on SUM of hits from last 3 attempts)
		set.seed(1110)#get it ?, miss after streak of 3
		r=c()
		for (k in 1:10000) r[k]=sim(100,p0=.376,p1=.45,p2=.508,p3=.742,type='moving.sum')$r
		mean(r)

     #Method=streak (hotness is base on streat of 0-3 in a row))
        set.seed(1110)
        r=c()
        for (k in 1:10000) r[k]=sim(100,p0=.376,p1=.45,p2=.508,p3=.742,type='streak')$r
        mean(r)
  ##############################################################################################################
  
 