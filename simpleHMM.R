#library(MASS)
#library(MCMCpack)
#############GENERATE SOME DATA
#Simulate an HMM
T=50#number of observations
K=2

lambda<-c(2,8)
G<-matrix(c(.8,.2,.2,.8),nrow = K,byrow = TRUE)
D<-solve(t(diag(K)-G+1),rep(1,K))
######
path<-sample(1:K,size = 1,prob = D)
for(i in 2:T){
        path[i]<-sample(1:K,size = 1,prob = G[path[i-1],])
}
#######
path.lambda<-sapply(path,function(x) lambda[x])
obs<-sapply(path.lambda,function(x) rpois(n = 1,lambda = x))
#######



#######---------------------------------------------------
#############################################################
#############################################################
########################################COMPUTE P############
cmpt.P<-function(OBSERVATIONS,LAMBDA){
        df<-data.frame(expand.grid(OBSERVATIONS,LAMBDA))
        df$Dpois<-apply(df,1,function(x)dpois(x[1],x[2]))
        B<-matrix(df$Dpois,ncol = length(LAMBDA),byrow = F)
        p.matrices<-lapply(data.frame(t(B)),diag)
        return(p.matrices)
}
#############################################################
#############################################################
########################################COMPUTE ALPHA########
#This function is depricated 
#cmpt.ALPHA<-function(MC.Dist,MC.tpm,P.matrices){
#        Y<-list()
#        Y[[1]]<-MC.Dist%*%P.matrices[[1]]
#        for(i in 2:T){
#                Y[[i]]<-Y[[i-1]]%*%MC.tpm%*%P.matrices[[i]]
#        }
#        alpha_next<-data.frame(
#                matrix(unlist(lapply(Y,function(x) x/sum(x))),#
#                       ncol = K,byrow = TRUE,
#                       dimnames = list(c(),dnames)))
#        return(alpha_next)
#}
cmpt.ALPHA<-function(MC.Dist,MC.tpm,P.matrices){
        w<-c(sum(MC.Dist%*%P.matrices[[1]]))
        phi<-matrix(0,ncol = K,nrow = T)
        phi[1,]<-(MC.Dist%*%P.matrices[[1]])/w[1]
        for(t in 2:T){
                w[t]<-w[t-1]*(phi[t-1,]%*%MC.tpm%*%P.matrices[[t]]%*%matrix(1,nrow = K))
                phi[t,]<-(w[t-1]/w[t])*phi[t-1,]%*%MC.tpm%*%P.matrices[[t]]
        }
                return(phi)
}




#############################################################table(factor(path,levels=c("c1","c2","c3")))
#############################################################
#####################SIMULATE PATH###########################
cmpt.PATH<-function(ALPHA,MC.tpm){
        sample_path=rep(0,times = T)
        sample_path[T]<-sample(1:K,1,prob = alpha[T,])
        for(i in (T-1):1){
                pr<-alpha[i,]*MC.tpm[,sample_path[i+1]]
                pr<-pr/sum(pr)
                sample_path[i]<-sample(1:K,1,prob = pr)
        }
        return(sample_path)
}
#############################################################
#############################################################
#############################################################
############################UPDATE T.P.M ####################
Update.TP<-function(NU,PATH){
        #tr.cnt counts each transaction from the PATH
        tr.cnt<-table(factor(PATH[1:(T-1)],levels = dnames),
                      factor(PATH[2:T],levels = dnames))
        nu_next<-NU+tr.cnt
        #tpm<-t(apply(nu_next,1,function(x) rdirichlet(n=1,x)))
        #colnames(tpm)<-rownames(tpm)<-dnames
        #return(list(tpm,nu_next))
        return(nu_next)
}
#############################################################
#######################REGIMES CONTRIBUTIONS#################
regime_cntrbutn<-function(PATH,OBSERVATIONS,TAU){
        regime<-data.frame(STATE=PATH,OBS=OBSERVATIONS)
        regime<-cbind(regime,
                        data.frame(matrix(0,ncol = K,nrow = T,
                                                  dimnames = list(c(),paste("tau",1:K,sep = "")))))
        for(i in 1:K)regime[,(i+2)]<-(regime$STATE>=i)*TAU[i]
        M<-t(apply(regime,1,function(x) c(x[2],x[3:(2+K)]/sum(x[3:(2+K)]))))
        N<-data.frame(t(apply(M,1,function(x) 
                rmultinom(n = 1,size = x[1],prob = x[2:(K+1)]))))
        colnames(N)<-paste("regime",1:K,sep = "")
        regime<-cbind(regime,N)
        return(regime)
}

#############################################################
#############################################################
##########################Update Emission shpae.rate#########
Update.EMS<-function(PATH,OBSERVATIONS,TAU){
        regimes_cnt<-regime_cntrbutn(PATH = PATH,OBSERVATIONS = OBSERVATIONS,TAU=TAU)
        sum.jt<-unname(apply(regimes_cnt[,(3+K):(2+(2*K))],2,sum))
        N.j<-   unname(apply(regimes_cnt[,3:(K+2)],2,function(x) sum(x>0)))
        c<-c(sum.jt,N.j)
        tau.inc<-t(matrix(c,nrow = 2,byrow = TRUE,dimnames = list(c("shape","rate"))))
        SH.R<-shape.rate+tau.inc
        return(SH.R)
}
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###################INITIALIZE THE PARAMETERS
###################Do it either manually or
###################Randomly


dnames<-paste("c",1:K,sep = "")
nu<-matrix(sample(10:100,size = K^2,replace = TRUE),nrow = K)
G<-t(apply(nu,1,function(e) rdirichlet(n = 1,alpha = e)))
colnames(G)<-rownames(G)<-dnames
#----
shape.rate<-data.frame(shape=sample(2:7,size = K,replace = TRUE),rate=sample(x = 1:2,size=K,replace = TRUE))
tau<-apply(X = shape.rate,MARGIN = 1,function(x) rgamma(n = 1,shape = x[1],rate = x[2]))
lambda<-cumsum(x = tau) #Emission Rate




###------------------------------------------------------------
iteration=100000
BURN_IN=1000
Glist<-list()
Llist<-list()

ptm <- proc.time()
for(i in 1:iteration){
        P<-cmpt.P(OBSERVATIONS = obs, LAMBDA = lambda)
        alpha<-cmpt.ALPHA(MC.Dist = D,MC.tpm = G,P.matrices= P)
        path<-cmpt.PATH(ALPHA = alpha,MC.tpm = G)
        #--
        nu<-Update.TP(NU = nu,PATH = path)
        G<-t(apply(nu,1,function(x) rdirichlet(n=1,x)))
        colnames(G)<-rownames(G)<-dnames
        #--
        shape.rate<-Update.EMS(PATH = path,OBSERVATIONS = obs,TAU = tau)
        tau<-apply(shape.rate,1,function(x) rgamma(1,shape = x[1],rate = x[2]))
        lambda<-cumsum(tau)
        #--
        D<-solve(t(diag(K)-G+1),rep(1,K))
        #D<-rdirichlet(n = 1,alpha =  as.vector(unname(prop.table(table(path)))))
        if(i%%1000==0)print(paste("@ iteration ",str(i)))
        if(i>BURN_IN){
                Glist[[i-BURN_IN]]<-G
                Llist[[i-BURN_IN]]<-lambda
        }
       
}
proc.time()-ptm




