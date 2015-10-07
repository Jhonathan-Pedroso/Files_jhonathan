sim.M.poli<- function (n.gen,n.ind,n.loc) {
  n.gam=2
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
  }
  return (M.f)
}

A.meu <- function (K,n,ploi) {
  
  p=colSums(K)/(ploi*nrow(K))
  Kp=(matrix(0,nrow(K),ncol(K)))
  
  for (i in 1:nrow(K)) {
    for (j in 1:ncol(K)){
      Kp[i,j]=(K[i,j]-n*(p[1,j]))/sqrt(n*p[1,j]*(1-p[1,j]))
    }
  }
  A<-tcrossprod(Kp)/ncol(K)
  A<-A+diag(nrow(A))*1e-04
}
A.vit= function (K,n) {
  freq=colSums(K)/(ploi*nrow(K))
  Fa=t(matrix(freq,ncol(K),nrow(K)))
  
  W=K-n*Fa
  cor=n*sum(freq*(1-freq))
  Av=W%*%t(W)/cor
  return (Av)
}
Wa.vit= function (K) {
  freq=colSums(K)/(ploi*nrow(K))
  Fa=t(matrix(freq,ncol(K),nrow(K)))
  
  W=K-5*Fa
  return (W)
}

require(rrBLUP)
library(ggplot2)
library(gridExtra)
library(agricolae)

ploi=4
A=A.vit(M,5)
W=Wa.vit(M)
n.ana=200
n.faixa.herd=4
herd=2
for (herd in 1:n.faixa.herd) {
  for (i in 1:n.ana) { 
    
    if (herd==1) {
      h2 <- 0.2  #heritability  
    }
    if (herd==2) {
      h2 <- 0.4  #heritability  
    }
    if (herd==3) {
      h2 <- 0.6  #heritability  
    }
    if (herd==4) {
      h2 <- 0.8  #heritability  
    }
    
    y <- 1000 +g + rnorm(nrow(A),mean=0,sd=sqrt((1-h2)/h2*var(g)))
    
    if (i==1) {
      ans <- mixed.solve(y,K=A)
      accuracy <- cor(g,ans$u)
      Vu = ans$Vu
      h=cov(ans$u,y)/var(y)  
      MSE.g=sum(c(ans$u-g)^2)/length(g)
      
      acc.vec=matrix(c(0),n.ana,1)
      Vu.vec=matrix(c(0),n.ana,1)
      h.vec=matrix(c(0),n.ana,1)
      MSE.g.vec=matrix(c(0),n.ana,1)
      acc.vec[i,1]=accuracy
      Vu.vec[i,1]=Vu
      h.vec[i,1]=h
      MSE.g.vec[i,1]=MSE.g
    }
    if (i!=1) {
      ans <- mixed.solve(y,K=A)
      accuracy <- cor(g,ans$u)
      Vu = ans$Vu
      h=cov(ans$u,y)/var(y)
      MSE.g=sum(c(ans$u-g)^2)/length(g)
      
      acc.vec[i,1]=accuracy
      Vu.vec[i,1]=Vu
      h.vec[i,1]=h
      MSE.g.vec[i,1]=MSE.g
    }
  }
  
  if (herd==1) {
    tiff(filename = "herd2.tiff",width=1100, height=700, units="px")
    m=qplot(c(h.vec), geom="histogram",xlim=c(0,1)) + theme_bw() 
    m + geom_histogram(aes(fill = ..count..)) +
      scale_fill_gradient("Contagem", low = "grey", high = "gold") +
      labs(x = "Frequência alélica", y = "") +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,2)")
    dev.off()
    
    x=c(ans$u);y=c(g);
    dat=data.frame(A=x,B=y)
    VGG2= ggplot(data=dat, aes(x=A,y=B)) 
    VGG2= VGG2 + theme_bw() + 
      labs(x = "VGG estimado", y = "VGG paramétrico") +
      geom_point(colour = "gold", size = 3.5,shape = 18) +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,2)") + stat_smooth(method="lm", se=FALSE,color="black")
    
    x=c(acc.vec);y=c(MSE.g.vec);
    dat=data.frame(A=x,B=y)
    MC2= ggplot(data=dat, aes(x=A,y=B)) 
    MC2= MC2 + theme_bw() + 
      labs(x = "Correlações", y = "MSE") + xlim(0,1) + ylim(0,max(y))+
      geom_point(colour = "gold", size = 3.5,shape = 18) +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,2)")
    
    Vu2=round(mean(c(Vu.vec)),2)
    
  }
  if (herd==2) {
    
    tiff(filename = "herd4.tiff",width=1100, height=700, units="px")
    m=qplot(c(h.vec), geom="histogram",xlim=c(0,1)) + theme_bw() 
    m + geom_histogram(aes(fill = ..count..)) +
      scale_fill_gradient("Contagem", low = "grey", high = "green") +
      labs(x = "Frequência alélica", y = "") +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,4)")
    dev.off()
    
    x=c(ans$u);y=c(g);
    dat=data.frame(A=x,B=y)
    VGG4= ggplot(data=dat, aes(x=A,y=B)) 
    VGG4= VGG4 + theme_bw() + 
      labs(x = "VGG estimado", y = "VGG paramétrico") +
      geom_point(colour = "green", size = 3.5,shape = 18) +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,4)") + stat_smooth(method="lm", se=FALSE,color="black")
    
    x=c(acc.vec);y=c(MSE.g.vec);
    dat=data.frame(A=x,B=y)
    MC4= ggplot(data=dat, aes(x=A,y=B)) 
    MC4= MC4 + theme_bw() + 
      labs(x = "Correlações", y = "MSE") + xlim(0,1) + ylim(0,max(y))+
      geom_point(colour = "green", size = 3.5,shape = 18) +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,4)")
    
    Vu4=round(mean(c(Vu.vec)),2)
    
  }
  if (herd==3) {
    tiff(filename = "herd6.tiff",width=1100, height=700, units="px")
    m=qplot(c(h.vec), geom="histogram",xlim=c(0,1)) + theme_bw() 
    m + geom_histogram(aes(fill = ..count..)) +
      scale_fill_gradient("Contagem", low = "grey", high = "blue") +
      labs(x = "Frequência alélica", y = "") +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,6)")
    dev.off()
    
    x=c(ans$u);y=c(g);
    dat=data.frame(A=x,B=y)
    VGG6= ggplot(data=dat, aes(x=A,y=B)) 
    VGG6= VGG6 + theme_bw() + 
      labs(x = "VGG estimado", y = "VGG paramétrico") +
      geom_point(colour = "blue", size = 3.5,shape = 18) +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,6)") + stat_smooth(method="lm", se=FALSE,color="black")
    
    x=c(acc.vec);y=c(MSE.g.vec);
    dat=data.frame(A=x,B=y)
    MC6= ggplot(data=dat, aes(x=A,y=B)) 
    MC6= MC6 + theme_bw() + 
      labs(x = "Correlações", y = "MSE") + xlim(0,1) + ylim(0,max(y))+
      geom_point(colour = "blue", size = 3.5,shape = 18) +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,6)")
    
    Vu6=round(mean(c(Vu.vec)),2)
    
  }
  if (herd==4) {
    tiff(filename = "herd8.tiff",width=1100, height=700, units="px")
    m=qplot(c(h.vec), geom="histogram",xlim=c(0,1)) + theme_bw() 
    m + geom_histogram(aes(fill = ..count..)) +
      scale_fill_gradient("Contagem", low = "grey", high = "red") +
      labs(x = "Frequência alélica", y = "") +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,8)")
    dev.off()
    
    x=c(ans$u);y=c(g);
    dat=data.frame(A=x,B=y)
    VGG8= ggplot(data=dat, aes(x=A,y=B)) 
    VGG8= VGG8 + theme_bw() + 
      labs(x = "VGG estimado", y = "VGG paramétrico") +
      geom_point(colour = "red", size = 3.5,shape = 18) +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,8)") + stat_smooth(method="lm", se=FALSE,color="black")
    
    x=c(acc.vec);y=c(MSE.g.vec);
    dat=data.frame(A=x,B=y)
    MC8= ggplot(data=dat, aes(x=A,y=B)) 
    MC8= MC8 + theme_bw() + 
      labs(x = "Correlações", y = "MSE") + xlim(0,1) + ylim(0,max(y))+
      geom_point(colour = "red", size = 3.5,shape = 18) +
      theme(text = element_text(size=20,face="bold",colour="grey20"),
            axis.text.x = element_text(colour="grey20",size=20,angle=360,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain")) +
      ggtitle("h² (0,8)")
    
    Vu8=round(mean(c(Vu.vec)),2)
  }
  
  if (herd==n.faixa.herd) {
    
    tiff(filename = "VGG.tiff",width=1100, height=700, units="px")
    
    grid.arrange(VGG2, VGG4, VGG6, VGG8, ncol=2)
    
    dev.off()
    
    tiff(filename = "MSE_cor.tiff",width=1100, height=700, units="px")
    
    grid.arrange(MC2, MC4, MC6, MC8, ncol=2)
    
    dev.off()
    
    Vu.table=matrix(c(Vu2,Vu4,Vu6,Vu8),byrow=T)
    rownames(Vu.table) <- c("Vu2","Vu4","Vu6","Vu8")
    save(Vu.table,file="Vu.table.RData")
  }
}
