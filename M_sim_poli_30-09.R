sim.M.poli<- function (n.gen,n.ind,n.loc) {
  n.gam=2
  salto=c(seq(from=1, to=2000, by=19))
  G1=as.matrix(c(1,0,0,0,0))  #AAAA
  G2=as.matrix(c(0,1,0,0,0))  #AAAa
  G3=as.matrix(c(0,0,1,0,0))  #AAaa
  G4=as.matrix(c(0,0,0,1,0))  #Aaaa
  G5=as.matrix(c(0,0,0,0,1))  #aaaa
  
  gam1=as.matrix(c(1,0,0))  #AA
  gam2=as.matrix(c(0,1,0))  #Aa
  gam3=as.matrix(c(0,0,1))  #aa
  for (loc in 1:n.loc) {
    for (gen in 1:n.gen) {
      if (gen==1) {
        P=rmultinom(n.ind,1,prob=c((1/5),(1/5),(1/5),(1/5),(1/5)))  ##It is the matrix of the individuals genotypes, where the columns represents de individuals  
      }
      if (gen!=1) {
        P=Prog
      }
      
      for (j in 1:n.ind) {
        index=cbind(as.matrix(P[,j])==G1,as.matrix(P[,j])==G2,as.matrix(P[,j])==G3,
                    as.matrix(P[,j])==G4,as.matrix(P[,j])==G5)   ##index is the logical matrix to find the individuals genotypes
        
        #Testing if all elements of a vector is true is done by "all()"
        if (j==1) {
          gam=matrix(0,3,n.ind*n.gam) #gam is the gametes, where the columns represent the individuals
        }
        #The "& i==" it is to make sure that the gametic will be given for the proper genotype; it is necessary to match precisely in the sequence
        
        for (i in 1:5) {
          if (all(index[,i]) & i==1) {   
            gam[,c(j,(n.ind+j))]=gam1
          }
          if (all(index[,i]) & i==2) {   
            gam[,c(j,(n.ind+j))]=rmultinom(n.gam,1,prob=c(0.5,0.5,0))
          } 
          if (all(index[,i]) & i==3) {   
            gam[,c(j,(n.ind+j))]=rmultinom(n.gam,1,prob=c((1/6),(4/6),(1/6)))
          }
          if (all(index[,i]) & i==4) {   
            gam[,c(j,(n.ind+j))]=rmultinom(n.gam,1,prob=c(0,0.5,0.5))
          }
          if (all(index[,i]) & i==5) {   
            gam[,c(j,(n.ind+j))]=gam3
          }
        }
      }
      
      ##Defining the progenies
      
      Prog=matrix(0,5,n.ind)
      vec=matrix(1:(n.gam*n.ind),n.gam,n.ind,byrow=F)
      for (i in 1:(n.ind)) {
        if (all(as.matrix(gam[,c(vec[1,i])])==gam1) & all(as.matrix(gam[,c(vec[2,i])]==gam1))) {
          Prog[,i]=G1
        }
        if (all(as.matrix(gam[,c(vec[1,i])])==gam1) & all(as.matrix(gam[,c(vec[2,i])]==gam2))) {
          Prog[,i]=G2
        }
        if (all(as.matrix(gam[,c(vec[1,i])])==gam1) & all(as.matrix(gam[,c(vec[2,i])]==gam3))) {
          Prog[,i]=G3
        }
        if (all(as.matrix(gam[,c(vec[1,i])])==gam2) & all(as.matrix(gam[,c(vec[2,i])]==gam1))) {
          Prog[,i]=G2
        }
        if (all(as.matrix(gam[,c(vec[1,i])])==gam2) & all(as.matrix(gam[,c(vec[2,i])]==gam2))) {
          Prog[,i]=G3
        }
        if (all(as.matrix(gam[,c(vec[1,i])])==gam2) & all(as.matrix(gam[,c(vec[2,i])]==gam3))) {
          Prog[,i]=G4
        }
        if (all(as.matrix(gam[,c(vec[1,i])])==gam3) & all(as.matrix(gam[,c(vec[2,i])]==gam1))) {
          Prog[,i]=G3
        }
        if (all(as.matrix(gam[,c(vec[1,i])])==gam3) & all(as.matrix(gam[,c(vec[2,i])]==gam2))) {
          Prog[,i]=G4
        }
        if (all(as.matrix(gam[,c(vec[1,i])])==gam3) & all(as.matrix(gam[,c(vec[2,i])]==gam3))) {
          Prog[,i]=G5
        }  
      }
      
      if (gen==n.gen) {
        
        M=matrix(0,n.ind,1)
        for (i in 1:n.ind) {
          logic <- all(Prog[,i]==G1)
          if (logic) {
            M[i,1] <-4
          }
          logic <- all(Prog[,i]==G2)
          if (logic) {
            M[i,1] <-3
          }
          logic <- all(Prog[,i]==G3)
          if (logic) {
            M[i,1] <-2
          }
          logic <- all(Prog[,i]==G4)
          if (logic) {
            M[i,1] <-1
          }
          logic <- all(Prog[,i]==G5)
          if (logic) {
            M[i,1] <-0
          }
        }
        if (loc==1) {
          M.f=M
        }
        if (loc!=1) {
          M.f=cbind(M.f,M)
        }
      }
    }
    print(c(loc))
    if (any(loc==salto)) {
      save(M.f,file="M_500g_8.RData")
    }
  }
  return (M.f)
}

setwd("/home/jhonathan/Results")
M=sim.M.poli(500,250,1000)