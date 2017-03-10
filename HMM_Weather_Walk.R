L       <- 3          # number of days
x       <- c(2,1,2)   # 2= bike, 1=walk. The length of this should be L.
trans   <- matrix(c(.75, .50, .25, .50),nr=2)           # trans(i,j):Prob of going from state i to j
emit    <- matrix(c(.70, .20, .30, .80),nr=2)           # emit[i,j]-Prob of emitting j (bike or walk) from state i
epsilon <- 1.0
B       <- 500        # burn-in
S       <- 1000       # number of samples
K       <- 10         # sample spacing
lik     <- vector(length=(B+S*K))                       # this will hold the likelihood for each iteration
sunny.count <- rep(0,L)                                 # sunny.count[i] = number of times day i was sunny
#    on the sampled iterations
# Initialize Z's randomly (why is this okay?)

z       <- rbinom(L, 1, .5)+1                           # 2=sun, 1=rain
##############################

# Iterate the updates
#=====================
for(t in 1:(B+S*K)) {
        # Update each day in turn
        # Start with the first day
        prob.sun      <- (.5*trans[2,z[2]]*emit[2,x[1]])/
                         (.5*trans[2,z[2]]*emit[2,x[1]]+.5*trans[1,z[2]]*emit[1,x[1]])
        u             <- runif(1)
        if(u<prob.sun)  {
                z[1] <- 2           # sunny
        }
        if(u>=prob.sun) {
                z[1] <- 1           # rainy
        }
        
        #   Update days in middle
        #-------------------------
        for(i in 2:(L-1)) {
                num      <- trans[z[i-1],2]*trans[2,z[i+1]]*emit[2,x[i]];
                denom    <- trans[z[i-1],1]*trans[1,z[i+1]]*emit[1,x[i]] + num;
                prob.sun <- num / denom;
                u=runif(1)
                if(u<prob.sun){
                        z[i] <- 2   # sunny
                }
                if(u>=prob.sun){
                        z[i] <- 1   # rainy
                }
        } # end iteration over days
        
        #End with last day
        #------------------
        
        num           <- trans[z[L-1],2]*epsilon*emit[2,x[L]]
        denom         <- trans[z[L-1],1]*epsilon*emit[1,x[L]] + num;
        prob.sun      <- num / denom;
        
        u             <-runif(1)
        if(u<prob.sun){
                z[L]        <- 2  # sunny
        }
        if(u>=prob.sun){
                z[L]        <- 1  # rainy
        }
        
        # Calculate likelihood
        # --------------------
        lik[t] <- epsilon
        #multiply by Pr(X_i|z_i)*Pr(z_i|z_i-1), for i from 2 to L
        for(i in 2:L) {
                lik[t] <- lik[t]*emit[z[i],x[i]]*trans[z[i-1],z[i]]
        }
        
        #multiply by Pr(X_1|z_1)*Pr(Z_1|Begin)
        lik[t] <- lik[t]*emit[z[1],x[1]]*.5
        
        if( t>B && t%%K==0){  #sample the iteration
                # Check whether each day is sunny
                sunny.count=sunny.count+z-1
        }
} # end "for each iteration t"

# Estimate probability that each day is sunny
# ===========================================
sunny.prob <- sunny.count/S
print("sunny.prob=")
print(sunny.prob)