#' Utinity Function Based on Expected Decision Strategy (E.st).
#'
#' This utinity function returns the probability to select one out of two-armed bandits with Bernoulli-outcomes (favorable and unfavorable) for the next (n+1)-th individual. 
#' E.st considerates the making-decision in a population with the homogeneously expected decision attitudes, which actually calculates the expected value of beta posterior distribution. 
#' @param x    the informed informaiton of n individuals (n=0,1,2,...) with consisted of 4 integers including counts of A$success, A$failure, B$success and B$failure. 
#' @details    A and B are examples of two treatment arms in the context, resulting into a binary outcomes success or failure after one patient to be treated. The four integeral outcomes are written as a vector.  
#' @details    In terms of every patient with respecting themseleves'decision attitudes, each decision process selecting either A or B treatment is a randomly probabilistical process which depends the output returned by their selected utinity function.   
#' @return prob_A   the probability to select A-arm, which is consisted of values 1, 0.5 or 0. 
#' @keywords NA
#' @export
#' @examples
#' AS<-13
#'# the successful counts of A treatment
#' AF<-5
#'#the failures of A treatment
#' BS<-2
#'#the successful counts of B treatment
#' BF<-1
#'#the failures of B treatment
#' N<- AS+AF+BS+BF
#'#the total number of patients treated (n)
#' E.st_utinity(x=c(AS,AF,BS,BF))
#'# Decision to how much probability to select arm A for the (n+1)-th patient with an expected attitude. 


E.st_utinity<- function(x){ 
  x. <- x+1
  a <- x.[1]/(x.[1]+x.[2])
  b <- x.[3]/(x.[3]+x.[4])
  return(prob_A=(sign(a-b)+1)/2)
}


#' Utinity Function Based on Targeting Decision Strategy (T.st).
#' 
#' T.st considerates the making-decision in a population with target (optimistic/pessimistic) attitudes, and calculates the probability of favorable rate more than a targeting value. 
#' The optimism hope their targeting values higher than the higher expected value of Beta posterior out of two. In contrast, the pessimism hope their targeting values less likely than the higher expected one.
#' @param x        the informed informaiton of n individuals (n=0,1,2,...) with consisted of 4 integers including A$success, A$failure, B$success and B$failure counts. 
#' @param w        a degree of attitude of an individual,and is to parameterize his or her targeting value.  
#' @param w > 0    on behalf of an optimistic individual.
#' @param w < 0    on behalf of a pessimistic individual.
#' @details  A and B are examples of two treatment arms in the context, resulting into a binary outcomes success or failure after one patient to be treated. The four integeral outcomes are written as a vector.  
#' @details  In terms of every patient with respecting themseleves'decision attitudes, each decision process select either A or B treatment is a randomly probabilistical process which depends the output returned by their selected utinity function.   
#' @details  w range from -1 to 1. 
#' @return prob_A   the probability to select arm-A, which is consisted of values 1, 0.5 or 0. 
#' @keywords NA
#' @export
#' @examples
#' 
#' AS<-13
#'# the successful counts of A treatment
#' AF<-5
#'#the failures of A treatment
#' BS<-2
#'#the successful counts of B treatment
#' BF<-1
#'#the failures of B treatment
#' N<- AS+AF+BS+BF
#'#the total number of patients treated (n)
#' w=0.5
#'# the degree of the (n+1)-th individual's attitude who is one of the optimism. 
#' T.st_utinity(x=c(AS,AF,BS,BF),w=0.5)
#'# Decide to how much probability to select arm A for the (n+1)-th patient with an optimistic attitude. 

T.st_utinity <- function(x,w){ 
  x. <- x+1
  a.exp <- x.[1]/(x.[1]+x.[2])
  b.exp <- x.[3]/(x.[3]+x.[4])
  tmp<-max(a.exp,b.exp)
  if(w>0){
    target<-tmp+(1-tmp)*w
  } else {
    target<-tmp*(1+w)
  }
  a <- pbeta(target,x.[1],x.[2],lower.tail=FALSE)
  b <- pbeta(target,x.[3],x.[4],lower.tail=FALSE)
  return(prob_A=(sign(a-b)+1)/2)
}



#' Enumeration of All Possible Four Integral Two-armed Bernoulli-outcomes. 
#'
#' This Function is to enumerate all possible two-armed Bernoulli-outcomes with the form of 2X2 tables When N patients to be treated. 
#' @param N        the number of patients to be treated.
#' @param X0       the initial information with setting c(0,0,0,0)
#' @return X=Xsm   a list of length N+1, each of which is a list of matrices. The values of matrices are 0 throughout.
#' @return M       the list of the output of 2-part of composition of intergers (n=0, 1,...,N) 
#' @return Z       a list with representing separately counts of A&success, A&failure, B&succsess, and B&failure; and the length of each list equals to the number of all possible combinations of Bernoulli-outcomes when N patients to be treated. 
#' @return sk      to get the same attributes of Xsm.  
#' @details   Decision process for the first patient starts with information 0, where no previous patients treated and so only one possible Bernoulli-outcoms with X0=c(0,0,0,0). 
#' @details   xsimple() function used here to enumerate all 2-part composition of N patients to be treated, on a simplex lattice (2,n=N). And this function is required to attach the package "combinat". 
#' @details   ZAS, ZAF,ZBS,and ZBF separately indicate the counts of  A&Success,A&Failure,B&Success, B&Failure.
#' @keywords NA
#' @export
#' @examples
#install.packages("combinat")
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' attributes(out)


serialTables <- function(N,X0=rep(0,4)){
  M <- list()
  for(i in 0:N){
    M[[i+1]] <- as.matrix(xsimplex(2,i)) 
  }
  Xsm <- ZAS <- ZAF <- ZBS <- ZBF <- list()
  ZAS[[1]] <- matrix(X0[1],1,1)
  ZAF[[1]] <- matrix(X0[2],1,1)
  ZBS[[1]] <- matrix(X0[3],1,1)
  ZBF[[1]] <- matrix(X0[4],1,1)
  Xsm[[1]] <- list(matrix(0,1,1))
  for(i in 1:N){
    n <- i+1
    Xsm[[n]] <- list()
    ZAS[[n]] <- ZAF[[n]] <- ZBS[[n]] <- ZBF[[n]] <- list()
    for(j in 1:n){
      Xsm[[n]][[j]] <- matrix(0,j,n-j+1)
      ZAS[[n]][[j]] <- matrix(rep((j-1):0,n-j+1),j,n-j+1)
      ZAF[[n]][[j]] <- (j-1) - ZAS[[n]][[j]]
    }
    for(j in 1:n){
      ZBS[[n]][[j]] <- t(ZAS[[n]][[n+1-j]])
      ZBF[[n]][[j]] <- t(ZAF[[n]][[n+1-j]])
    }
  }
  Z <- list(unlist(ZAS),unlist(ZAF),unlist(ZBS),unlist(ZBF)) 
  sk <- attr(unlist(as.relistable(Xsm),recursive=TRUE),"skeleton")
  return(list(X=Xsm,M=M,Z=Z,sk=sk))
}


#' Calculate Exact Probability of a 2X2 Table with Bernoulli-outcomes
#'
#' This function is to calculate the exact probability of every possible 2X2 table with two-armed Bernoulli-outcomes enumerated in the serialTables () when N patients to be treated.
#' @param s1     a list of probabilities to select A-arm for the all possible 2X2 tables based on E.st_utinity () or T.st_utinity ().
#' @param p      a vector of true success rates of two arms, a and b. 
#' @return       a list of exact probabilities per 2X2 table. 
#' @keywords NA
#' @export
#' @examples
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' ABSF<-sapply(out$Z, unlist)
#'
#' #Probability to select the A-arm based on utinity function of E.st: E.st_utinity()
#' E.st_prob.A<-apply(ABSF,1,E.st_utinity)
#' E.table.probability<-table.prob(relist( E.st_prob.A,out$sk),p=c(0.8,0.6))
#' E.table.probability.v<-unlist(E.table.probability)
#' E.table.probability.v
#'
#' #Probability to select the A-arm based on utinity function of T.st: T.st_utinity() 
#' T.st_prob.A<-apply(ABSF,1,T.st_utinity,w=0.5)
#' T.table.probability<-table.prob(relist(T.st_prob.A,out$sk),p=c(0.8,0.6))
#' T.table.probability.v<-unlist(T.table.probability)
#' T.table.probability.v
table.prob <- function(s1,p){
  p.A <- p[1]
  p.B <- p[2]
  ret <- s1
  ret[[1]][[1]] <- 1
  N <- length(s1)-1
  for(i in 1:N){
    n <- i+1
    ret[[n]] <- lapply(ret[[n]],"*",0)
    for(j in 1:i){
      tmp.select <- s1[[i]][[j]]
      dims <- dim(tmp.select)
      # (1) A & Success
      J <- j+1;xspan <- 1:dims[1];yspan <- 1:dims[2];
      ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + tmp.select * ret[[i]][[j]] * p.A
      # (2) A & Failure
      J <- j+1;xspan <- 2:(dims[1]+1);yspan <- 1:dims[2];
      ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + tmp.select * ret[[i]][[j]] * (1-p.A)
      # (3) B & Success
      J <- j;xspan <- 1:dims[1];yspan <- 1:dims[2];
      ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + (1-tmp.select) * ret[[i]][[j]]* p.B
      # (4) B & Failure
      J <- j;xspan <- 1:dims[1];yspan <- 2:(dims[2]+1);
      ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + (1-tmp.select) * ret[[i]][[j]] * (1-p.B)
    }
  }
  ret
}



#' Average Fraction of Favorable Outcomes of Two-armed Bernoulli-outcomes.
#' 
#' Calculate the average success fraction of two arms for every 2X2 table numerated in the serialTables () when N patients to be treated. 
#' @param x      a vector of 4 elements inlcuding A$Success, A$Failure, B$Success, and B$Failure.
#' @details   A and B are examples of two arms, which resulting in favorable outcomes (success) and unfavorable outcomes (failures).
#' @return       average of success fraction
#' @keywords NA
#' @export
#' @examples
#' AS<-13
#' #  the successful counts of A treatment
#' AF<-5
#' #  the failures of A treatment
#' BS<-2
#' #  the successful counts of B treatment
#' BF<-1
#' #  the failures of B treatment
#' N<- AS+AF+BS+BF
#' success.frac(x=c(AS,AF,BS,BF))
#'
#'
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' ABSF<-sapply(out$Z, unlist)
#' success.fraction<-apply(ABSF,1,success.frac)
#' #  the average fraction of successes per table.
#' success.fraction

success.frac <- function(x){
  if(sum(x)==0){
    0.5
  }else{
    fl <- x[2]+x[4]
    suc <- x[1]+x[3]
    suc/(suc+fl)
  }
  
}



#' Calculate the Fraction of A-arm Selection of Two-armed Bernoulli-outcomes.
#' 
#' Calculate the fraction of  A_armed for every 2X2 table numerated in the serialTables () when N patients to be treated. 
#' @param x     a vector of 4 elements including A$Success, A$Failure, B$Success, and B$Failure.
#' @details  A and B are examples of two arms, which resulting in favorable outcomes (success) and unfavorable outcomes (failure).
#' @return      the fraction of A-arm selected
#' @keywords NA
#' @export
#' @examples
#' AS<-13
#' #  the successful counts of A treatment
#' AF<-5
#' #  the failures of A treatment
#' BS<-2
#' # the successful counts of B treatment
#' BF<-1
#' #the failures of B treatment
#' N<- AS+AF+BS+BF
#' A.frac(x=c(AS,AF,BS,BF))
#' 
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' ABSF<-sapply(out$Z, unlist)
#' A.fraction<-apply(ABSF,1,A.frac)
#' # the fraction of A-arm selected per table.
#' A.fraction

A.frac <- function(x){
  if(sum(x)==0){
    0.5 
  }else{
    A.cont <- x[1]+x[2]
    Total.cont<- sum(x)
    A.cont/Total.cont
  }
}



#' Weighed Average Favorable Outcomes Rate per Individual
#' 
#' Weighted mean success rate per individual, named as overall success rate, since it calculates the sum of the multiplication of the probability and average success rate for all possible 2X2 tables when N patients. 

#' @param xv      a vector of average success fraction of two arms in all possible tables.
#' @param pv      a vector of occurence probabilities of all tables
#' @param sk      a list of list structure of matrices indicating all possible Bernoulli-outcomes of two arms.
#' @return        a vector of overall success rate of a series of patients n=0,1,2,...,N.
#' @keywords NA
#' @export
#' @examples
 #'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' ABSF<-sapply(out$Z, unlist)
#'# all possible 2X2 tables
#'
#' E.st_prob.A<-apply(ABSF,1,E.st_utinity)
#' #Probability to select the A-arm based on utinity function of E.st.
#' E.table.probability<-table.prob(relist(E.st_prob.A,out$sk),p=c(0.8,0.6))
#' E.table.probability.v<-unlist(E.table.probability)
#' success.fraction<-apply(ABSF,1,success.frac)
#' success.fraction[1]<-0.5
#' OSR(xv=success.fraction,pv=E.table.probability.v,sk=out$sk)
#'
#' # Probability to select the A-arm based on utinity function of T.st.
#' T.st_prob.A<-apply(ABSF,1,T.st_utinity,w=0.5)
#' T.table.probability<-table.prob(relist(T.st_prob.A,out$sk),p=c(0.8,0.6))
#' T.table.probability.v<-unlist(T.table.probability)
#' success.fraction<-apply(ABSF,1,success.frac)
#' success.fraction[1]<-0.5
#' OSR(xv=success.fraction,pv=T.table.probability.v,sk=out$sk) 
OSR <- function(xv,pv,sk){
  wv <- xv * pv
  wl <- relist(wv,sk)
  sapply(wl,function(x)sum(unlist(x)))
}


#' Overall Fraction of A-arm Selection  per Individual
#' 
#' Average the fraction of A-arm selection per individual.

#' @param fv     a vector of the A-arm selection fraction out of two-armed Bernoulli in all possible tables.
#' @param pv     a vector of occurence probabilities of all tables
#' @param sk     a list of list structure of matrices with indicating all possible Bernoulli-outcomes when N patients treated.
#' @return       a vector of mean of A-arm selection fraction of a series of patients n=0,1,2,...,N.
#' @keywords NA
#' @export
#' @examples
#'library(combinat) 
#' N <- 10
#' out <- serialTables(N)
#' ABSF<-sapply(out$Z, unlist)
#' #all possible 2X2 tables
#'
#' E.st_prob.A<-apply(ABSF,1,E.st_utinity)
#' #Probability to select the A-arm based on utinity function of E.st.
#' E.table.probability<-table.prob(relist(E.st_prob.A,out$sk),p=c(0.8,0.6))
#' E.table.probability.v<-unlist(E.table.probability)
#' A.fraction<-apply(ABSF,1,A.frac)
#' A.fraction[1]<-0.5
#' O_A_frac(fv=A.fraction,pv=E.table.probability.v,sk=out$sk)
#'
#' #Probability to select the A-arm based on utinity function of T.st.
#' T.st_prob.A<-apply(ABSF,1,T.st_utinity,w=0.5)
#' T.table.probability<-table.prob(relist(T.st_prob.A,out$sk),p=c(0.8,0.6))
#' T.table.probability.v<-unlist(T.table.probability)
#' A.fraction<-apply(ABSF,1,A.frac)
#' A.fraction[1]<-0.5
#' O_A_frac(fv=A.fraction,pv=T.table.probability.v,sk=out$sk) 
O_A_frac<- function(fv,pv,sk){
  wv <- fv * pv
  wl <- relist(wv,sk)
  sapply(wl,function(x)sum(unlist(x)))
}


   
   #' Graph of OSRs and A-arm Selection Fraction Based on homogeneous E.st and T.st.
   #' 
   #' Visulize the OSRs and A-arm selection Fraction of homogeneous E.st and homogeneous T.st when N patients treated in a plot. 
   #' @param  d       data of matrix with two colums, separately store the values based on the E.st and T.st
   #' @param  title   the name of the plot. 
  #' @return          the 2-D plot  
   #' @keywords NA
   #' @export
   #' @examples
   #' library(combinat)
   #' N <- 20
   #' out<-serialTables(N)
   #' ABSF<-sapply(out$Z, unlist)  
   #'  
   #' A.fraction<-apply(ABSF,1,A.frac)
   #' A.fraction[1] <- 0.5
   #' Success.fraction <-success.fraction<-apply(ABSF,1,success.frac)
   #' Success.fraction[1] <- 0.5
   #'
   #' E.st_prob.A<-apply(ABSF,1,E.st_utinity)
   #' #Probability to select the A-arm based on utinity function of E.st: E.st_utinity()
   #' a<-0.8
   #' #ture success rate of A-arm treatment
   #' b<-0.6
   #' #ture success rate of B-arm treatment
   #' E.table.probability<-table.prob(relist(E.st_prob.A,out$sk),p=c(a,b))
   #' #exact probability per 2X2 table based on utinity function of E.st: E.st_utinity()
   #' 
   #' w=0.5
   #' #homogeneous optimistic attitude with w=0.5
   #' T.st_prob.A<-apply(ABSF,1,T.st_utinity,w=0.5)
   #' #Probability to select the A-arm based on utinity function of T.st: T.st_utinity()
   #' a<-0.8
   #' #ture success rate of A-arm treatment
   #' b<-0.6
   #' #ture success rate of B-arm treatment
   #' T.table.probability<-table.prob(relist(T.st_prob.A,out$sk),p=c(a,b))
   #' #exact probability per 2X2 table based on utinity function of T.st: T.st_utinity()
   #' 
   #'#Comparison of A-arm Fraction and OSRs of homogeneous E.st and T.st
   #' E.table.prob.v<- unlist(E.table.probability)
   #' E.mean.selectionA.per.n<-O_A_frac(fv=A.fraction,pv=E.table.prob.v,sk=out$sk) 
   #' E.mean.success.per.n<-OSR(xv=success.fraction,pv=E.table.prob.v,sk=out$sk) 
   #' 
   #' T.table.prob.v<- unlist(T.table.probability)
   #' T.mean.selectionA.per.n<-O_A_frac(fv=A.fraction,pv=T.table.prob.v,sk=out$sk) 
   #' T.mean.success.per.n<-OSR(xv=success.fraction,pv=T.table.prob.v,sk=out$sk) 
   #'
   #' par(mar=c(4,4,3,1),mfrow=c(1,2))
   #' Vis_OSR_A.fra(d=cbind(E.mean.selectionA.per.n,T.mean.selectionA.per.n),title="A-arm Fraction of E.st vs T.st")
   #' Vis_OSR_A.fra(d=cbind(E.mean.success.per.n,T.mean.success.per.n),title="OSRs of E.st vs T.st")
 
Vis_OSR_A.fra<-function(d,title){
  matplot(d,type="b",pch=19,xlab=" ",ylab=" ",font.axis=4,cex.axis=1.2)
  mtext(title,side=3,cex=1.8,font=4,line=0.5)
  mtext(title,side=1,cex=1.6,font=4,line=2.5)
  mtext("Rate",side=2,cex=1.6,font=4,line=2.5)
  legend("bottomright",c("E.st","T.st"),col=1:2,text.col=1:2,cex=1.2,text.font=4,box.lwd=2,pch=19)
}
      
   
   
#' 2-D Grid of Three Variables
#' 
#' Make the 2-D matrix to store the probabilities of all 2X2 tables when N patients treated, columns and rows are A-arm selection fraction and sucess rate separately. . 
#' @param  S      the same structure with Xsm,storing the number of favorable outcomes at all possible 2X2 tables when a series of N patients to be treated.
#' @param  A      the same structure with Xsm,storing the number of A-arm treatment at all possible 2X2 tables when a series of N patients to be treated
#' @param  p      the same structure with Xsm,storing the probabilities of all possible 2X2 tables when a series of N patients to be treated.
#' @details   S,A,p are special list structure, each of which represents one variable based on all possible combinations of 2X2 tables when N patients treated. And the length of list equals to N+1.  
#' @return        a list of matrices consisted of probabilities with rows corresponding to the probabilities vary along with A-arm fraction selection and columns corresponding to the probabilities vary along with overall success rates.   
#' @keywords NA
#' @export
#' @examples
#' library(combinat)
#' N <- 20
#' out<-serialTables(N)
#' ABSF<-sapply(out$Z, unlist)  
#' #'  
#' A.cnt <- relist((ABSF[,1]+ABSF[,2]),out$sk)
#' A.fraction<-apply(ABSF,1,A.frac)
#' A.fraction[1] <- 0.5
#' Success.cnt <- relist(ABSF[,1]+ABSF[,3],out$sk)
#' Success.fraction <-success.fraction<-apply(ABSF,1,success.frac)
#' Success.fraction[1] <- 0.5
#'
   
#' E.st_prob.A<-apply(ABSF,1,E.st_utinity)
#' #Probability to select the A-arm based on utinity function of E.st: E.st_utinity()
#' a<-0.8
#' #ture success rate of A-arm treatment
#' b<-0.6
#' #ture success rate of B-arm treatment
#' E.table.probability<-table.prob(relist(E.st_prob.A,out$sk),p=c(a,b))
#' #exact probability per 2X2 table based on utinity function of E.st: E.st_utinity()
#' grid.exp.ls<-dist.2D(Success.cnt,A.cnt,E.table.probability)
#' str(grid.exp.ls)
#' 
#' w=0.5
#' #homogeneous optimistic attitude with w=0.5
#' T.st_prob.A<-apply(ABSF,1,T.st_utinity,w=0.5)
#' #Probability to select the A-arm based on utinity function of T.st: T.st_utinity()
#' a<-0.8
#' #ture success rate of A-arm treatment
#' b<-0.6
#' #ture success rate of B-arm treatment
#' T.table.probability<-table.prob(relist(T.st_prob.A,out$sk),p=c(a,b))
#' #exact probability per 2X2 table based on utinity function of T.st: T.st_utinity()
#' grid.target.ls<-dist.2D(Success.cnt,A.cnt,T.table.probability) 
#' str(grid.target.ls)


   dist.2D <- function(S,A,p){
     N <- length(p)
     grid.prob.series <- list()
     for(i in 1:N){
       grid.prob.series[[i]] <- matrix(0,i,i)
       p.v <- unlist(p[[i]])
       s.v <- unlist(S[[i]])
       a.v <- unlist(A[[i]])
       L <- length(p.v)
       for(j in 1:L){
         grid.prob.series[[i]][s.v[j]+1,a.v[j]+1] <- grid.prob.series[[i]][s.v[j]+1,a.v[j]+1] + p.v[j]
       }
     }
     return(grid.prob.series)
   }
   
#' Visulization of Probability, A-arm Selection Fraction and Overall Success Rate.
#' 
#' Visulization probabilities against the A-arm fraction and Overall success rate plane via persp()
#' @param g       grid matrices outputed from the function dist.2D ()
#' @param title   the name of the strategy. E.st or T.st  
#' @param zlims   the minimum and maximum z-axis.
#' @details x-axis indicates A-arm fraction and y-axis indicates the overall success rate, and z-axis indicates the probability. 
#' @return  perspective plot with surface     
#' @keywords NA
#' @export
#' @examples
#' library(combinat)
#' N <- 50
#' out<-serialTables(N)
#' ABSF<-sapply(out$Z, unlist)  
#'  
#' A.cnt <- relist((ABSF[,1]+ABSF[,2]),out$sk)
#' A.fraction<-apply(ABSF,1,A.frac)
#' A.fraction[1] <- 0.5
#' Success.cnt <- relist(ABSF[,1]+ABSF[,3],out$sk)
#' Success.fraction <-success.fraction<-apply(ABSF,1,success.frac)
#' Success.fraction[1] <- 0.5
#'
   
#' E.st_prob.A<-apply(ABSF,1,E.st_utinity)
#' #Probability to select the A-arm based on utinity function of E.st: E.st_utinity()
#' a<-0.8
#' #ture success rate of A-arm treatment
#' b<-0.6
#' #ture success rate of B-arm treatment
#' E.table.probability<-table.prob(relist(E.st_prob.A,out$sk),p=c(a,b))
#' #exact probability per 2X2 table based on utinity function of E.st: E.st_utinity()
#' grid.exp.ls<-dist.2D(Success.cnt,A.cnt,E.table.probability)
#' 
#' w=0.5
#' #homogeneous optimistic attitude with w=0.5
#' T.st_prob.A<-apply(ABSF,1,T.st_utinity,w=0.5)
#' #Probability to select the A-arm based on utinity function of T.st: T.st_utinity()
#' a<-0.8
#' #ture success rate of A-arm treatment
#' b<-0.6
#' #ture success rate of B-arm treatment
#' T.table.probability<-table.prob(relist(T.st_prob.A,out$sk),p=c(a,b))
#' #exact probability per 2X2 table based on utinity function of T.st: T.st_utinity()
#' grid.target.ls<-dist.2D(Success.cnt,A.cnt,T.table.probability) 
#'  
#' zlims<-range(c(grid.exp.ls[[N+1]]),c(grid.target.ls[[N+1]]))
#' par(bg="white",cex=1,col="black",mar=c(2,3.5,3,1))
#' Vis_3d(g=grid.exp.ls[[N+1]],title="Homogeneous E.st",zlims=zlims)
#' Vis_3d(g=grid.target.ls[[N+1]],title="Homogeneous T.st",zlims=zlims)
   Vis_3d<-function(g,title,zlims){
     persp(g,r=10, theta =200, phi = 20, expand = 0.7, col= "springgreen" ,
           xlab="Success Rate",ylab="A-arm Fraction",zlab="Probability", ticktype="detailed",font.axis=4,cex.axis=1.1,
           main=title,cex.main=1.8,cex.lab=1.5,font.main=4,font.lab=4,border=NULL,zlim=c(zlims))
   }
   
   
   
#' Rescale the Colors Based on Raw Dataset. 
#' 
#' Color-scale based on raw data points with your visualized colors. 
#' @param obs            a verctor of raw data   
#' @param neg.cols       the colors to be assigned for negatives.        
#' @param pos.cols       the colors to be assigned for positives.   
#' @details rgb() is used for rescaling the colors and "gdata" package is necessary to be installed. 
#' @return  the re-scaled colors corresponding the raw data points.   
#' @keywords NA
#' @export
#' @examples
#' load("mean.success.per.n.exp.100.201.RData")
#' load("mean.success.per.n.wpbeta.100.201.RData")  
#' N <- 100
#' ps<-combn(seq(0,1,by=0.01),2)
#' w<-seq(from=-1,to=1,by=0.01)
#' Difs<-array(0,c(length(w),length(ps[1,]),N))
#'    for(i in 1:N){
#'      Difs[,,i]<-t(t(mean.success.per.n.wpbeta[,,i+1])-mean.success.per.n.exp[i+1,])
#'    }
#'    w.id<-c(21,61,141,181)
#'    #the positition of w values in the vector w, such as w=0.8,-0.4,0.4,or 0.8
#'    N.id<-c(10,30,100)
#'    Difs.<-Difs[w.id,,N.id]
#'    cols<-color_scale(obs=c(Difs.),neg.cols=c("blue","white"),pos.cols=c("white","red"))
   #'    par(mfrow=c(1,1),mar=c(3,3,3,3))
   #'    image(matrix(seq(0,1,length.out=length(c(Difs.))),ncol=1),col=cols[order(c(Difs.))], xlab=" ", ylab=" ",xaxt="n",yaxt="n")
   #'    axis(3,at=c(0.3),labels=c("Negatives"),cex.axis=1.4,font=4,las=1,lwd=4,tick=FALSE,line=-0.2)
   #'    axis(3,at=c(0.85),labels=c("Positives"),cex.axis=1.4,font=4,las=1,lwd=4,tick=FALSE,line=-0.2)
   #'    se<-range(c(Difs.))
   #'    axis(3,at=c(0,0.5,1),labels=round(c(min(Difs.),median(Difs.),max(Difs.)),3),cex.axis=1.1,font=4,las=1,line=-0.2,lwd=4,col.ticks="green",col="green")
   #'    colbar<-seq(0,1,length.out=length(cols))
   #'    axis(3,at=colbar[length(which(Difs.<0))],labels=c("0"),cex.axis=1.4,font=4,las=1,lwd=4,line=-0.2,col.ticks="green")
   
   color_scale<-function(obs,neg.cols,pos.cols){
     dif.val<-obs/(max(obs)-min(obs))
     neg.id<-which(dif.val<0)
     pos.id<-which(dif.val>0)
     zero.id<-which(dif.val==0)
     dif.val.neg<-dif.val[neg.id]-min(obs)/(max(obs)-min(obs))
     dif.val.pos<-dif.val[pos.id]-min(obs)/(max(obs)-min(obs))
     dif.val.zero<-dif.val[neg.id]
     
     col_fun.neg <- colorRamp(neg.cols)
     dif.val.neg.1<-seq(0,1,length.out=length(dif.val.neg))
     rgb_cols.neg <- col_fun.neg(c(dif.val.neg.1))  
     cols.neg<- rgb(rgb_cols.neg, maxColorValue = 256) 
     cols.neg.1<-dif.val.neg
     cols.neg.1[order(cols.neg.1)]<-cols.neg
     
     
     col_fun.pos <- colorRamp(pos.cols)
     dif.val.pos.1<-seq(0,1,length.out=length(dif.val.pos))
     rgb_cols.pos <- col_fun.pos(c(dif.val.pos.1))  
     cols.pos<- rgb(rgb_cols.pos, maxColorValue = 256) 
     cols.pos.1<-dif.val.pos
     cols.pos.1[order(cols.pos.1)]<-cols.pos
     
     
     
     cols<-rep("black",length(dif.val))
     cols[neg.id]<-cols.neg.1 
     cols[pos.id]<-cols.pos.1
     return(cols)
   }
   
   
   
     
   
#' Visulization of the Benefits of Homogeneous T.st with Few Typical w Cases and E.st.
#' 
#' Visulization of difference of OSRs based on homogeneous T.st and E.st with few cases of w values, and on the condition of multiple ture success pairs, a, b. 
#' @param D          a matrix of difference of OSRs between homogeneous T.st and E.st with rows (the length of w), columns (a,b  pairs). 
#' @param interval   the interval of sequence of success rate a/b range from [0,1].
   #' @param N.p        the number of patient treated.
   #' @param cols.mat    the same structure with D with keeping the corresponding color coded. 
#' @details Actually we've already calculated the OSRs based on the homogeneous T.st and E.st under the conditions w = {-0.99,-0.98,...,0,0.01,...,0.98,0.99}, N= 1,2,.,100, for (a,b) = {(pa,pb)|pa,pb = {0.01,0.02,.,0.98,0.99},a>=b}. 
#' @details Each process of calculation is the same as the examples in the function OSR () with given N, one pair of (a, b), and w. 
#' @details Differences with fixed N and w on the 5050 (a,b) pairs are coded colored in the triangle (gdata package used here), the blue indicates negatives and the red inicates positives.
#' @details We saved the simulated data based on E.st and T.st with data form .RData separately in this package, and we can load two datasets to make the image.
#' @return  trangle-contour image coded  with blue-red color. 
#' @keywords NA
#' @export
#' @examples
#' load("mean.success.per.n.exp.100.201.RData")
#' load("mean.success.per.n.wpbeta.100.201.RData")
#' library(gdata)
#' N <- 100
#' interval=0.01
#' ps<-combn(seq(0,1,by=interval),2)
#' w<-seq(from=-1,to=1,by=0.01)
#' Difs<-array(0,c(length(w),length(ps[1,]),N))
#'    for(i in 1:N){
#'      Difs[,,i]<-t(t(mean.success.per.n.wpbeta[,,i+1])-mean.success.per.n.exp[i+1,])
#'    }
#'    w.id<-c(21,61,141,181)
#'    N.id<-c(10,30,100)
#'    Difs.<-Difs[w.id,,N.id]
#'    cols<-color_scale(obs=c(Difs.),neg.cols=c("blue","white"),pos.cols=c("white","red"))
#' cols.<-array(cols,c(dim(Difs.)))
#' par(mfrow=c(length(N.id),length(w.id)),lwd=2)
#' for(i in seq_along(N.id)){
  #'    N.p<-N.id[i]
  #'    D<-Difs.[,,i]
  #'    cols.mat<-cols.[,,i]
  #'    Vis_image.hom(D=D,cols.mat=cols.mat,N.p=N.p,interval=interval)
  #' }
  

   Vis_image.hom<-function(D,cols.mat,N.p,interval){
    
      for(j in seq_along(w.id)){
        lengs<-length(seq(0,1,by=interval))
           if(j==1){
             par(mar=c(0.1,2.2,2.2,0.1))
             mat<-D[j,]
             ps.pairs.mean.dif<-matrix(NA,lengs,lengs)
             lowerTriangle(ps.pairs.mean.dif)<-mat
             ps.pairs.mean.dif[101,]<-NA
             ps.pairs.mean.dif[,1]<-NA
             image(ps.pairs.mean.dif, xlab=" ", ylab=" ", xaxt="n",yaxt="n",col=cols.mat[j,][order(mat)])
             axis(2,at=c(0,1),labels=c("0","1"),cex.axis=2,font=4,las=3,line=-0.9,tick=FALSE)
             axis(3,at=c(0,1),labels=c("0","1"),cex.axis=2,font=4,las=1,line=-0.9,tick=FALSE)
             mtext(paste0("w = ",w[w.id[j]]),side=3,cex=2,font=4,line=0.15)
             mtext(paste0("N = ",N.p),side=2,cex=2,font=4,line=0.15)
             arrows(0.04,0.97,0.95,0.97,length=0.05)
             arrows(0.02,0.04,0.02,0.97,length=0.05)
             mtext("(a)",side=3,cex=2,font=4,line=-2.5)
             mtext("(b)",side=2,cex=2,font=4,line=-2.5)
             contour(ps.pairs.mean.dif,add=TRUE,labcex=1,drawlabels=FALSE,col="green")
           }else{
             par(mar=c(0.1,2.2,2.2,0.1))
             mat<-D[j,]
             ps.pairs.mean.dif<-matrix(NA,lengs,lengs)
             lowerTriangle(ps.pairs.mean.dif)<-mat
             ps.pairs.mean.dif[101,]<-NA
             ps.pairs.mean.dif[,1]<-NA
             image(ps.pairs.mean.dif, xlab=" ", ylab=" ", xaxt="n",yaxt="n",col=cols.mat[j,][order(mat)])
             mtext(paste0("w = ",w[w.id[j]]),side=3,cex=2,font=4,line=0.15)
             contour(ps.pairs.mean.dif,add=TRUE,labcex=1,drawlabels=FALSE,col="green")
           }
           
         }
       }


   
   
      
   

   #' Nine Classification in Comparison of OSRs of Different Types of Decision Attitudes Popualtion. 
   #' 
   #' Rgeardless of decision attitudes (w values), we builded nine classifications based on comparison of OSRs between homogeneous T.st(optimistic/pessimistic attitudes) and E.st population(expected attitudes). 
#' @param D.mat        a matrix of difference of OSRs between homogeneous T.st and E.st with rows (the length of w), columns (a,b pairs).  
#' @param w            a vector of the degree of homogeneoudecision attitudes. 
#' @details Nine classifications based on 2 w levels including -1< w <0 and 0<w<1, and 3 levels including OSRs of T.st - E.st  are all  non-negtives (+/0) or all negatives (-) with regardless of w values, or sometimes positives at particular w cases or sometimes negatives (+/-).
#' @details Nine possible classifications numerized as integers from 1 to 9.
#' @details  1 states that regardless of w, (optimistic/pessimistic) T.st population perform no worser than E.st
#' @details  2 states that regardless of w, the pessimist perform better or equal than E.st , however the optimist with outperforming E.st depend on w value.
#' @details  3 states that regardless of w, values, the pessimist perform best, but the optimist performs worest. 
#' @details  4 states that regardless of w, the optimist perform best, however E.st or the pessimist either better or worser depend on the w values.
#' @details  5 states that which type of population outperform depend on w values.
#' @details  6 states that the optimist are worest, and who, the E.st or the pessimist,perform better depend on w values. 
#' @details  7 states that regardless of w, the optimist perform no worser than E.st, and E.st perform better than the pessimist.
#' @details  8 states that regardless of w, the pessimist perform no worest, others depend on w values. 
#' @details  9 states that regardless of w, the E.st population perform T.st whatever optimistic or pessimistic attitudes 
#' @details  Actually, in terms of without caring w values, 1,4,7 classifications together indicate the optimist outperform E.st, 1,2,3 together indicate the pessimist outperform E.st, and 9 indicates the E.st population outperform better than T.st.  
#' @details  The loaded dataset in the examples is simulated based on the conditions that w = {-0.99,-0.98,...,0,0.01,...,0.98,0.99}, N= 1,2,.,100, for (a,b) = {(a,b)|a,b \in {0.01,0.02,.,0.98,0.99},a>=b}. 
#' @return  a matrix with the 5050 rows (number of a,b pairs) and 3 columns with separately representing a, b, and numeric index of classification (1,2,...,9).
#' @keywords NA
#' @export
#' @examples
#' load("mean.success.per.n.exp.100.201.RData")
#' load("mean.success.per.n.wpbeta.100.201.RData")
#' 
#' N <- 100
#' ps<-combn(seq(0,1,by=0.01),2)
#' w<-seq(from=-1,to=1,by=0.01)
#' Difs<-array(0,c(length(w),length(ps[1,]),N))
#'    for(i in 1:N){
#'      Difs[,,i]<-t(t(mean.success.per.n.wpbeta[,,i+1])-mean.success.per.n.exp[i+1,])
#'    }
#'   
#' index.100<-cls.fun(D.mat=Difs[,,N],w=w)
#' # When patients N=100
#' table(index.100[,3])  

   
   
  cls.fun<-function(D.mat,w){
     w.non_postive<-which(w>-1&w<0)
     w.positive<-which(w<1&w>0)
     Dif.1<-D.mat[w.non_postive,]
     Dif.2<-D.mat[w.positive,]
     type.id<-matrix(0,length(ps[1,]),3)
     for(i in seq_along(ps[1,])){
       type.id[i,1]<-ps[1,i]
       type.id[i,2]<-ps[2,i]
       if(all(round(Dif.1[,i],15)>=0)){
         if(all(round(Dif.2[,i],15)<0)){
           type.id[i,3]<-3
         }else if(all(round(Dif.2[,i],15)>=0)){
           type.id[i,3]<-1
         }else{
           type.id[i,3]<-2
         }
       }else if(all(round(Dif.1[,i],15)<0)){
         if(all(round(Dif.2[,i],15)<0)){
           type.id[i,3]<-9
         }else if(all(round(Dif.2[,i],15)>=0)){
           type.id[i,3]<-7 
         }else{
           type.id[i,3]<-8 
         }
       }else{
         if(all(round(Dif.2[,i],15)<0)){
           type.id[i,3]<-6
         }else if(all(round(Dif.2[,i],15)>=0))
         {
           type.id[i,3]<-4 
         }else{
           type.id[i,3]<-5
         }
       }
     } 
     return(type.id)
   }
 
   
   
   
   
  
#' Visulization of Benefits of Homogeneous T.st and E.st Regardless of w
#' 
#' Visulization of 9 classification after obtained numerical outputs from the above cls-fun(), and coded them with different colors in the triangle image. 
#' @param  index     a matrices with 5050 raws (number of a,b  pairs) and 3 columns (a,b,and numeric classification integers, 1,2,...,9).
#' @param  N.p         patients number N.
#' @param  interval   the interval of sequence of success rate a/b range from [0,1].
#' @details  The image is made in the example of this function is based on the simulated data ".RData" on the conditions that w = {-0.9995,-0.9990,...,0,0.0005,...,0.9990,0.9995}, N= 1,2,.,100, for (a,b) = {(a,b)|a,b = {0.01,0.02,.,0.98,0.99},a>=b}. 
#' @details  Differences of OSRs were calulated at each condition a,b pair (5050 pairs) and each w values (4001 values). And this dataset ".RData" loaded had been classified via function cls_fun(). 
#' @details "gdata" package used here, and color the blue indicates negatives and the red inicates positives in the trangle conditions of a, b pairs. 
#' @details We saved the simulated data based on E.st and T.st with data form .RData separately in this package, and we can load two datasets to make the image.
#' @return  trangle image with contour coded  with blue-red color. 
#' @keywords NA
#' @export
#' @examples
#' load("type.index.out.0.0005.RData", ex <- new.env())
#' names(ex)
#' N <- 100
#' interval<-0.01
#' ps<-combn(seq(0,1,by=interval),2)
#' w<-seq(-1,1,by=0.0005)
#' load("type.index.out.0.0005.RData")
#' Index.out<-type.index.out
#' 
#' table(Index.out[[N]][,3])
#' # 1    4    5    6    9 
#' # 1  828 3831  209  181 
#' #the table of 9 classification when patient N=100
#' b<-Index.out[[N]][which(Index.out[[N]][,3]==1),1]
#' b
#' a<-Index.out[[N]][which(Index.out[[N]][,3]==1),2]
#' a
#' # when true success rate a=1, and b=0, the classification 1 exists. 4 with 828 cases of a.b pairs  indicates the optimist outperformed E.st and 
#' # 9 with 181 pairs indicates E.st outperfomed T.st
#' 
#' library(gdata)
#' #visulized the above classifications with red, blue and gray separately.
#' N.id<-c(10,30,50,100)
#' par(mfrow=c(1,length(N.id)))
#' for(i in c(N.id)){
#' index<-Index.out[[i]]
#' N.p=i
#' Vis_image.cls(N.p=N.p,index=index,interval=interval)
#' }
#' 

  Vis_image.cls<-function(N.p,index,interval){
    par(mar=c(2.2,2.4,2.5,1))
    par(lwd=3) 
    type.index<-index[,3]
    
    #a=1 and b=0 is not included here
  type.index[type.index==4]<-1
  type.index[type.index==5]<-2
  type.index[type.index==6]<-2
  type.index[type.index==9]<-3
    
    lengs<-length(seq(0,1,by=interval))
    ps.types<-matrix(NA,lengs,lengs)
    lowerTriangle(ps.types)<-type.index
    ps.types[101,]<-NA
    ps.types[,1]<-NA
    
    image(ps.types, xlab=" ",ylab=" ", xaxt="n",yaxt="n",col=c("red","gray50","blue"))
    axis(1,at=c(0,0.25,0.5,0.75,1),labels=c("0","0.25","0.5","0.75","1"),cex.axis=1.8,font=4,las=1,lwd=1.5,line=-0.1)#,tick=FALSE)
    axis(2,at=c(0,0.25,0.5,0.75,1),labels=c("0","0.25","0.5","0.75","1"),cex.axis=1.8,font=4,las=3,lwd=1.5,line=-0.1)
    mtext("(a)",side=3,cex=1.8,font=4,line=-2.7)
    mtext("(b)",side=2,cex=1.8,font=4,line=-2.7)
    mtext(paste0("N= ", N.p),side=3,cex=2,font=4,line=0.3)
    arrows(0.04,0.97,0.95,0.97,length=0.05)
    arrows(0.02,0.04,0.02,0.97,length=0.05)
  }
 

   
   #' Distribution of Heterogeneity of Decision Attitudes
   #' 
   #' Density function of decision attitudes with heterogeneous w is assumed from the modified beta distribution.
   #' @param u,v     non-negative parameters of the Beta distribution.
   #' @details  u=v indicates the symmetric shape of w distribution, indicating equal chance pessimism and optimism in a population. 
   #' @details  The central zero indicates the homogeneity of expected attitude population. More flat the shape is, and more variance of heterogeneity of decision attitude is. 
   #' @return   the attitude degree w randomly generated from the distribution.
   #' @keywords NA
   #' @export
   #' @examples
   #' u=5
   #' v=5
   #' n=10000
   #' w_obs<-W(n,u,v)
   #' par(mfrow=c(1,1))
   #' hist(w_obs)
   W<-function(n,u,v){
     w<--1+2*rbeta(n,u,v)
     return(w)
   }
   
   #' The Model of Heterogeneity Decision Attitudes
   #' 
   #' Return the overall success rate (OSR) based on the model of heterogeneity decision attitudes (optimistic/pessimistic).
   #' @param  iter      the times of randomly generated from W distribution. 
   #' @param  N         the number of patients to be treated
   #' @param  a,b       the true success rates of two arms A and B treatment
   #' @param  u,v       two shape parameters of the distributions of decision attitudes w.  
   #' @return            the mean of iter times of OSR
   #' @keywords NA
   #' @export
   #' @examples 
   #'library(combinat) 
   #'iter<-50
   #' N<-10
   #' a<-0.8
   #' b<-0.6
   #' u<-5
   #' v<-5
   #' Het_func(iter,N,a,b,u,v)
   Het_func<-function(iter,N,a,b,u,v){
     out<-serialTables(N)
     AS <- out$Z[[1]]
     AF <- out$Z[[2]]
     BS <- out$Z[[3]]
     BF <- out$Z[[4]]
     ABSF <- cbind(AS,AF,BS,BF)
     FL <- ABSF[,2]+ABSF[,4]
     SUC <- ABSF[,1]+ABSF[,3]
     Success.frac <- SUC/(SUC+FL)
     Success.frac[1]<- 0.5
     legs<-sapply(lapply(out$sk,function(x)sapply(x,length)),sum)
     legs.cum<-cumsum(legs)
     
     mean.success.per.n.het<-matrix(0,iter,(N+1))
     prob_A<-list()
     prob_A[[1]] <- 0.5
     for(i in 1:iter){
       for(j in 1:N){
         #Number of tables corresponds to each individual
         mat<-ABSF[(legs.cum[j]+1):legs.cum[j+1],]
         w<-W(1,u=u,v=v)
         #w<--1+2*rbeta(1,u,v)
         mat. <- mat+1
         a.exp <- mat.[,1]/(mat.[,1]+mat.[,2])
         b.exp <- mat.[,3]/(mat.[,3]+mat.[,4])
         tmp<-pmax(a.exp,b.exp)
         if(w>0){
           target<-tmp+(1-tmp)*w
         } else {
           target<-tmp*(1+w)
         }
         a <- pbeta(target,mat.[,1],mat.[,2],lower.tail=FALSE)
         b <- pbeta(target,mat.[,3],mat.[,4],lower.tail=FALSE)
         prob_A[[j+1]] <- (sign(a-b)+1)/2
       }             
       ####calculate OSR
       table.prob.v<-table.prob(relist(unlist(prob_A),out$sk),p=c(a,b))
       table.prob.het.v<-unlist(table.prob.v)
       mean.success.per.n.het[i,]<-OSR(xv=Success.frac,pv=table.prob.het.v,sk=out$sk)
       
     }
     
     het.mean.osr<-apply(mean.success.per.n.het,2,mean)
     return(het.mean.osr)
   }
   
   

   
   #' Visulization of Benefits of heterogeneity T.st and E.st
   #' 
   #' Visulization of difference of OSRs based on heterogeneous T.st and E.st on the condition of multiple ture success pairs, a, b.  
   #' @param  het.D      a vector of the difference of OSRs between heterogenous T.st and E.st. 
   #' @param  U.v         the value of parameter of W distribution 
   #' @param  interval    the interval of sequence of success rate a/b range from [0,1].
   #' @param  col.v      a vector of the color elements coded re-scaled by raw data points via color_scale().
   #' @details  ".RData" on the conditions that,N= 1,2,.,100,for (a,b) = {(a,b)|a,b={0.05,0.1,.,0.90,0.95},a>=b}. 
   #' @details  w generated from distribution function W() with various parameters u=v=5, 30, 70 or 461.
   #' @return   The trangle image with colored contour.           
   #' @keywords NA
   #' @export
   #' @examples 
   #' N<-50
   #' interval<-0.05
   #' ps<-combn(seq(0,1,by=interval),2)
   #' u=v=c(5,30,70,461)
   #' load("mean.success.per.n.exp.100.ps210.RData")
   #' exp.50<-mean.success.per.n.exp[N+1,]
   #' load("het.OSRs.RData")
   #' str(het.OSRs)
   #' het<-sapply(het.OSRs,function(x)x[,N+1])
   #' het.dif<-het-exp.50
  
    #' cols<-color_scale(c(het.dif),neg.cols=c("blue","white"),pos.cols=c("white","red"))
   #' par(mfrow=c(1,1),mar=c(3,3,3,3))
   #' image(matrix(seq(0,1,length.out=length(c(het.dif))),ncol=1),col=cols[order(c(het.dif))], xlab=" ", ylab=" ",xaxt="n",yaxt="n")
   #' axis(3,at=c(0.2),labels=c("Negatives"),cex.axis=1.4,font=4,las=1,lwd=4,tick=FALSE,line=-0.2)
   #' axis(3,at=c(0.75),labels=c("Positives"),cex.axis=1.4,font=4,las=1,lwd=4,tick=FALSE,line=-0.2)
   #' se<-range(c(het.dif))
   #' axis(3,at=c(0,0.5,1),labels=round(c(min(het.dif),median(het.dif),max(het.dif)),3),cex.axis=1.1,font=4,las=1,line=-0.2,lwd=4,col.ticks="green",col="green")
   #' colbar<-seq(0,1,length.out=length(cols))
   #' axis(3,at=colbar[length(which(het.dif<0))],labels=c("0"),cex.axis=1.4,font=4,las=1,lwd=4,line=-0.2,col.ticks="green")
   #' cols.mat<-matrix(cols, ncol=ncol(het.dif), nrow=nrow(het.dif))
   
   #' library(gdata)
   #' par(mfrow=c(1,length(u)))
   #' for(i in seq_along(u)){
   #'   het.D<-het.dif[,i]
   #'   U.v<-u[i]
   #'   col.v<-cols.mat[,i]
   #'   Vis_image.het(het.D=het.D,U.v=U.v, interval=interval,col.v=col.v)
   #' }

   Vis_image.het<-function(het.D,U.v,interval,col.v) {
     par(mar=c(2.2,2.4,2.5,3))
     mat<-het.D
     lengs<-length(seq(0,1,by=interval))
     ps.pairs.mean.dif<-matrix(NA,lengs,lengs)
     lowerTriangle(ps.pairs.mean.dif)<-mat
     ps.pairs.mean.dif[21,]<-NA
     ps.pairs.mean.dif[,1]<-NA
     image(ps.pairs.mean.dif,xlab=" ", ylab=" ", xaxt="n",yaxt="n",col=col.v[order(mat)])
     axis(1,at=c(0,0.25,0.5,0.75,1),labels=c("0","0.25","0.5","0.75","1"),cex.axis=1.8,font=4,las=1,lwd=1.5,line=-0.1)#,tick=FALSE)
     axis(2,at=c(0,0.25,0.5,0.75,1),labels=c("0","0.25","0.5","0.75","1"),cex.axis=1.8,font=4,las=3,lwd=1.5,line=-0.1)
     mtext(paste0("u=v= ",U.v),side=3,cex=2,font=4,line=0.15)
     arrows(0.04,0.97,0.95,0.97,length=0.05)
     arrows(0.02,0.04,0.02,0.97,length=0.05)
     mtext(at=0.7,"(a)",side=3,cex=1.8,font=4,line=-2.7)
     mtext(at=0.3,"(b)",side=2,cex=1.8,font=4,line=-2.7)
     contour(ps.pairs.mean.dif,add=TRUE,labcex=1.2,drawlabels=FALSE,col="green")
   }  
       
      
      