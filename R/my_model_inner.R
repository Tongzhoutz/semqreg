  
    rm(list=ls())
  library(pacman)
  library(data.table)
  library(tidyverse)
  library(haven)
  library(mpoly)
  library(corpcor)
  library(purrr)
  library(ggplot2)
  library(quantreg)
  library(magrittr)
  library(rio)
  library(progress)
  library(CVXR)
  library(Rcpp)
  library(RcppEigen)
  library(RcppArmadillo)
  library(RcppNumerical)
  
  df <- import("first_stage_regs.dta",setclass = "data.table")[, .(person,year,utoty,uc,ua,age)]
  
  DT <- copy(df)  ## ordered by year, NOT person!
  
  DT <- DT[order(year,person)]   ### match data in matlab file
  
  T = 6
  N = nrow(DT)/T
  
  source("_load_Hermite.R")
  
  ggplot(DT[year==98,],aes(age)) + geom_histogram(binwidth=1) + theme_economist()+
    labs(x = "age_98",title="Age distribution of households in 1998")
  
  DT[, str_c("Herm.age_",1:5) := hermPolys_4( (.SD-DT[,mean(age)] )/DT[,sd(age)] ), .SDcol = c("age")  ]
  
  coeff_DT <- DT[ , .(coef  = coef(lm(utoty ~ Herm.age_1+Herm.age_2+Herm.age_3+Herm.age_4+Herm.age_5 -1))  )]
  kappa_DT <- DT[, .(kappa = predict(lm(utoty ~ Herm.age_1+Herm.age_2+Herm.age_3+Herm.age_4+Herm.age_5 -1)))]
  
  DT[, predict_c := predict( lm( uc ~ Herm.age_1 + Herm.age_2+Herm.age_3+Herm.age_4+Herm.age_5-1 )) ]
  
  DT[, uc_resid := uc - predict_c]
  DT <- data.table(DT ,kappa_DT)
  uc_resid = DT[["uc_resid"]]
  ### residuals of OLS of utoty and uc on (age1, age^2, age^3, age^4)
  DT[ , utoty_nor := utoty - kappa  ][, -"kappa"]  
  DT[,  uc_nor := uc - predict_c]
  ### Hermite polynomials of utoty_nor
  DT[, str_c("Herm.year.",1:4) := hermPolys_3( (.SD-DT[,mean(utoty_nor)] )/DT[,sd(utoty)] ), .SDcol = c("utoty_nor")  ]
  utoty_nor = DT[["utoty_nor"]]
  ### Some statistics later will be used
  meanY = DT[, mean(utoty_nor)]
  stdY = DT[,sd(utoty_nor)]
  meanAGE = DT[,mean(age)]
  stdAGE = DT[,sd(age)]
  
  ################################################
                       ### M-Hastings Algorithm ##
  ################################################
  
  ## parameters 
  
    ## orders of Hermite polynomials
  set.seed(1003)
  
  ## the following adds tuning parameters of model
  
  source("_tuning.R")
  
  ############################################################
  ######### initial conditions: U_t Given  U_{t-1} and AGE_t  
  ############################################################
  ## k1 = 3
  source("_Ut.R")
  
  ############################################################
  ######### initial conditions: U_1 Given age_1  
  ############################################################
   ##   k2 =3
  source("_U1.R")
  
  ############################################################
  ######### initial conditions: V_t Given V_{t-1} + AGE_t  
  ############################################################
  ## k3 = 2
  
  source("_Vt.R")
  
  ############################################################
  ######### initial conditions: V_1 Given age_1  
  ############################################################
  ##   k4 =2
  source("_V1.R")
  
  ############### simulation ################################
  
  Resqnew_Ut <- array(NA, c(nrow = (k1+1)*(k_year+1),Ntau,maxiter))
  Resqnew_U1 <- array(NA, c(nrow =  k2+1  , Ntau, maxiter ))
  Resqnew_Vt <- array(NA, c(nrow =  2*(k3+1), Ntau, maxiter))
  Resqnew_V1 <- array(NA, c(nrow =  k4+1, Ntau, maxiter))
  
  init <- matrix( rnorm(N*T) , nrow = N, ncol=T )
  Herm_age <- DT[, .(Herm.age_1,Herm.age_2,Herm.age_3,Herm.age_4,Herm.age_5)] %>% as.matrix
  #### Tail parameters
  
  b1_ut = 10
  bL_ut = 10
  b1_u1 = 10
  bL_u1 = 10
  b1_vt = 10
  bL_vt = 10
  b1_v1 = 10
  bL_v1 = 10
  
  
  ### Mat used to store densities of U_t
  Mat_1 = matrix(NA,N,1)
  Mat_15 <- matrix(NA,N,5)
  Mat_26 <- matrix(NA,N,5)
  
  source("post_mcmc.R")
    sourceCpp("my_model_cpp.cpp")
  init = matrix( rnorm(N*T),N,T )
  
  Obj_chain = matrix(NA, nrow = N, ncol = draws)
  Obj_chain[, 1] = post_mcmc(init)
  
   posterior_cpp2(init,utoty_nor,uc_resid,Herm_age,regressor_V1,regressor_U1, Mat_init_Vt ,Mat_init_V1 ,Mat_init_Ut,Mat_init_U1,b1_ut, bL_ut, b1_u1,bL_u1, b1_vt,bL_vt,b1_v1,bL_v1 ,a,b,tauList,meanY,stdY )
  ones <- matrix(1 , N, draws)
  ones_vec <- rep(1,N)
  
  Nu_chain1 <- ones * (   init[,1] %*% t( rep(1,draws)) )
  Nu_chain2 <- ones * (   init[,2] %*% t( rep(1,draws)) )
  Nu_chain3 <- ones * (   init[,3] %*% t( rep(1,draws)) )
  Nu_chain4 <- ones * (   init[,4] %*% t( rep(1,draws)) )
  Nu_chain5 <- ones * (   init[,5] %*% t( rep(1,draws)) )
  Nu_chain6 <- ones * (   init[,6] %*% t( rep(1,draws)) )
  acc1 <- matrix(0,N,draws)
  acc2 <- matrix(0,N,draws)
  acc3 <- matrix(0,N,draws)
  acc4 <- matrix(0,N,draws)
  acc5 <- matrix(0,N,draws)
  acc6 <- matrix(0,N,draws)
  acceptrate1 <- rep(NA,draws)
  acceptrate2 <- rep(NA,draws)
  acceptrate3 <- rep(NA,draws)
  acceptrate4 <- rep(NA,draws)
  acceptrate5 <- rep(NA,draws)
  acceptrate6 <- rep(NA,draws)
  last_ut <- matrix(NA,N,T)
  mat_b <- matrix(NA, nrow = maxiter, ncol = 8)
  mat_lik <- rep(NA,maxiter)  
  MatU= matrix(0, N*(T-1),12)
  MatU2 = matrix(0,Ntau,12)
  newObj <- rep(NA,N)
  r <- rep(NA,N)
  prob <- rep(NA,N)
  resid_e1 = matrix(NA,N,11)
  resid_eps = matrix(NA,N*T,11)
  resid_tot = matrix(NA, N*(T-1),11)
  ut_reg_Mstep = matrix(NA, nrow = N*(T-1),ncol = 12)
  vt_reg_Mstep = matrix(NA,nrow= N*(T-1),ncol=6)
  
  #rand1 = rnorm(N)
  #rand2 = rnorm(N)
  #rand3 = rnorm(N)
  #rand4 = rnorm(N)
  #rand5 = rnorm(N)### E-step 
  #rand6 = rnorm(N)#source("_main.R")
  #probl = matrix(runif(N*T),N,T)
  #main_model(a,b,1 )
  
  #obj_tune <- map2_dfc(tune_ab[,1],tune_ab[,2], ~ main_model(.x,.y,200)    )
  
    ptm <- proc.time()
  for ( iter in 1:maxiter){
    print(iter)
    
      obj = innerLoop( draws,N,T,tauList,var_prop1 , var_prop2 ,var_prop3 , var_prop4 ,
               var_prop5 , var_prop6 , Nu_chain1 , Nu_chain2 , Nu_chain3,
               Nu_chain4, Nu_chain5, Nu_chain6 , Obj_chain , acc1,acc2,
               acc3, acc4,acc5,acc6,utoty_nor, uc_resid, Herm_age, regressor_V1 ,
               regressor_U1 , Mat_init_Vt , Mat_init_V1, Mat_init_Ut, Mat_init_U1 ,
               b1_ut, bL_ut, b1_u1, bL_u1, b1_vt, bL_vt, b1_v1, bL_v1, a, b,c,d,e,
               meanY,stdY)  
    
      acceptrate = obj[[4]]
      print( c( mean(acceptrate[,1]),mean(acceptrate[,2]),  mean(acceptrate[,3]), mean(acceptrate[,4]), 
                mean(acceptrate[,5]), mean(acceptrate[,6]) ) )
   last_ut = obj[[1]] 
    
    ##### M-step
    ut_2T <- as.vector( last_ut[,2:T] )
    ut_lag <- as.vector(last_ut[,1:5])
    
    Herm_ut_regressor <- sapply(  (ut_lag - meanY)/stdY , hermPolys_3 ) %>% t
    colnames(Herm_ut_regressor) <- str_c("V",11:14)
    
    for (i in 1:4){
      for (j in 1:3){
        ut_reg_Mstep[,3*(i-1) + j ] = Herm_ut_regressor[, i ] * Herm_age_2T[,j]
      }
    }
    #DT_ut_M <- data.table(Herm_ut_regressor,Herm_age_2T)
    #ut_reg_Mstep <-  DT_ut_M[, .(
    #                  C1 = V11  * V21,
    #                  C2 = V11 * V22,
    #                  C3 = V11 * V23,
    #                  C4 = V11 * V24,
    #                  C5 = V12  * V21,
    #                  C6 = V12 * V22,
    #                  C7 = V12 * V23,
    #                  C8 = V12 * V24,
    #                  C9 = V13  * V21,
    #                  C10 = V13 * V22,
    #                  C11 = V13 * V23,
    #                  C12 = V13 * V24,
    #                  C13 = V14  * V21,
    #                  C14 = V14 * V22,
    #                  C15 = V14 * V23,
    #                  C16 = V14 * V24
    #                  )] %>% as.matrix
                      
      fit_ut <- rq( ut_2T ~ ut_reg_Mstep-1, tau = tauList, method = "fn"  )
      Resqnew_Ut[,,iter] <- coef(fit_ut)
      
      fit_u1 <- rq( last_ut[,1] ~ regressor_U1-1, tau = tauList, method = "fn" )
      Resqnew_U1[,,iter] <- coef(fit_u1)  
      
      #DT[, ut_stack := as.vector(last_ut)]
      #vt_2T <- DT[, ut_stack := as.vector(last_ut) ][year != 98, .(vt_2T = utoty_nor - ut_stack)] %>% 
        #       as.matrix
      vt_2T =  utoty_nor[(N+1):(N*T)]- ut_2T
      DT_vt_15 = utoty_nor[1:(N*(T-1))] - ut_lag
      #DT_vt_15 <- DT[year != 108, .( vt_15 = utoty_nor - ut_stack  )] %>% as.matrix
      DT_vt_15_reg <- sapply( (DT_vt_15-meanY)/stdY, hermPolys_2  ) %>% t
      #DT_vt_age <- DT[year != 98, .( Herm.age_1, Herm.age_2 )] %>% as.matrix
      for (i in 1:3){
      for (j in 1:2){
        vt_reg_Mstep[,2*(i-1) + j ] = DT_vt_15_reg[, i ] * Herm_age_2T[,j]
      }
    }
    #
      #vt_reg_Mstep <- data.table( 
      #               C1= DT_vt_15_reg[,1] * DT_vt_age[,1],
      #               C2= DT_vt_15_reg[,1] * DT_vt_age[,2],
      #               C3= DT_vt_15_reg[,2] * DT_vt_age[,1],
      #               C4= DT_vt_15_reg[,2] * DT_vt_age[,2],
      #               C5= DT_vt_15_reg[,3] * DT_vt_age[,1],
      #               C6= DT_vt_15_reg[,3] * DT_vt_age[,2]) %>% as.matrix
        
      fit_vt <- rq(vt_2T ~ vt_reg_Mstep-1, tau = tauList, method = "fn")          
      Resqnew_Vt[,,iter] <- coef(fit_vt)       
            
      v1_dep = utoty_nor[1:N] - last_ut[,1]
      fit_v1 <- rq(v1_dep ~ regressor_V1-1, tau = tauList, method = "fn")
      Resqnew_V1[,,iter] <- coef(fit_v1)
      
      #### Normalization
      
      ##back_ut = mean(Resqnew_Ut[4,,iter])
     ## back_vt = mean(Resqnew_Vt[3,,iter])
      Resqnew_Ut[,,iter] <- Resqnew_Ut[,,iter] - apply( Resqnew_Ut[,,iter],1,mean )  ## demean by 3 rows
     Resqnew_Ut[4,,iter] <- Resqnew_Ut[4,,iter]+ 0.4 
      Resqnew_Vt[,,iter] <- Resqnew_Vt[,,iter] - apply( Resqnew_Vt[,,iter],1,mean )  ## demean by 3 rows
      
      #Resqnew_Vt[3,,iter] <- Resqnew_Vt[3,,iter] + 1.0 
      Resqnew_Ut[1,,iter] <- Resqnew_Ut[1,,iter] - ( (1 - tauList[11])/bL_ut - tauList[1]/b1_ut )*matrix(1,1,11)
   
      Resqnew_Vt[1,,iter] <- Resqnew_Vt[1,,iter] - ( (1 - tauList[11])/bL_vt - tauList[1]/b1_vt )*matrix(1,1,11)
      ### tail parameters
      Vect_b1_ut <-   ut_2T -   ut_reg_Mstep %*% Resqnew_Ut[,1,iter]
      Vect_bL_ut <-   ut_2T -   ut_reg_Mstep %*% Resqnew_Ut[,Ntau,iter]
      b1_ut <-  - sum( Vect_b1_ut <= 0) / sum( Vect_b1_ut *(Vect_b1_ut <=0) )
      bL_ut <- sum( Vect_bL_ut >= 0) / sum( Vect_bL_ut * (Vect_bL_ut >=0) )
      
      Vect_b1_u1 <-   last_ut[,1] - regressor_U1 %*% Resqnew_U1[,1,iter]
      Vect_bL_u1 <-   last_ut[,1] - regressor_U1 %*% Resqnew_U1[,Ntau,iter]
      b1_u1 <-  - sum( Vect_b1_u1 <= 0) / sum( Vect_b1_u1 *(Vect_b1_u1 <=0) )
      bL_u1 <- sum( Vect_bL_u1 >= 0) / sum( Vect_bL_u1 * (Vect_bL_u1 >=0) )
      
      Vect_b1_vt <-  vt_2T - vt_reg_Mstep%*% Resqnew_Vt[,1,iter]
      Vect_bL_vt <-  vt_2T - vt_reg_Mstep %*% Resqnew_Vt[,Ntau,iter]
      b1_vt <-  - sum( Vect_b1_vt <= 0) / sum( Vect_b1_vt *(Vect_b1_vt <=0) )
      bL_vt <- sum( Vect_bL_vt >= 0) / sum( Vect_bL_vt * (Vect_bL_vt >=0) )
      
      Vect_b1_v1 <-   v1_dep - regressor_V1 %*% Resqnew_V1[,1,iter]
      Vect_bL_v1 <-   v1_dep - regressor_V1 %*% Resqnew_V1[,Ntau,iter]
      b1_v1 <-  - sum( Vect_b1_v1 <= 0) / sum( Vect_b1_v1 *(Vect_b1_v1 <=0) )
      bL_v1 <- sum( Vect_bL_v1 >= 0) / sum( Vect_bL_v1 * (Vect_bL_v1 >=0) )
      
      ### for the next outer loop 
      
          ## for V_t 
           Mat_init_Vt <- Resqnew_Vt[,,iter]
          
          ## for V_1
           Mat_init_V1 <- Resqnew_V1[,,iter]
          
          ## for U_t
           Mat_init_Ut <- Resqnew_Ut[,,iter]
           
           ## for U_1
           Mat_init_U1 <- Resqnew_U1[,,iter]
           
           
           ## tail parameters
           mat_b[iter,1] <- b1_ut
           mat_b[iter,2] <- bL_ut
           mat_b[iter,3] <- b1_u1
           mat_b[iter,4] <- bL_u1
           mat_b[iter,5] <- b1_vt
           mat_b[iter,6] <- bL_vt
           mat_b[iter,7] <- b1_v1
           mat_b[iter,8] <- bL_v1
           
           ## persistence 
        #  for (i in 1:3){
        #    for (j in 0:2 ){
        #           MatU[,3*(i-1) + j + 4] = (i*Herm_ut_regressor[,i]/stdY)*Herm_age_2T[,(j+1)]
        #    }
        #  }
           
         #  Vect = quantile(ut_lag, tauList)
        #   U_prepare = hermPolys_2( ( as.vector(Vect)- meanY)/stdY )/stdY
        #   age_prepare = hermPolys_2( (meanAGE- meanAGE)/stdAGE)
         #  for (i in 1:3){
        #     for (j in 0:2){
        #       MatU2[,(3*i + j+1)] = i*U_prepare[,i] * age_prepare[j+1]
        #     }
        #   }
        #  print(MatU2 %*% Mat_init_Ut)
           ## mean log-likelihood
        
         #  print(apply(MatU %*% Mat_init_Ut,2,mean))
           mat_lik[iter] <-  mean(log(posterior_cpp2(last_ut,utoty_nor,uc_resid,Herm_age,regressor_V1,regressor_U1, Mat_init_Vt ,Mat_init_V1 ,Mat_init_Ut,Mat_init_U1,b1_ut, bL_ut, b1_u1,bL_u1, b1_vt,bL_vt,b1_v1,bL_v1 ,a,b,c,d,e,tauList,meanY,stdY )))
           if (is.na(mat_lik[iter] == TRUE )){
             stop('NA happens')}
           print(mat_lik[iter])
           
           #print( sd(last_ut[,6]))
           Obj_chain = obj[[3]]
           
   Nu_chain1 <- obj[[2]][,,1] 
   Nu_chain2 <- obj[[2]][,,2]
   Nu_chain3 <-  obj[[2]][,,3]
   Nu_chain4 <-  obj[[2]][,,4]
   Nu_chain5 <-  obj[[2]][,,5]
   Nu_chain6 <-  obj[[2]][,,6]
   
           ## 
  }  
  
    proc.time()-ptm
  mean( mat_lik[(maxiter/2) : maxiter])
maxiter = 500
mat_lik = A1_c[[3]]
df_gg = data.frame( x = 1:maxiter, y = mat_lik  )
  df_ggg = ggplot(df_gg, aes(x = x, y = mat_lik )) +
    geom_line(size=1) +
    labs(x = "iterations", y="likelihood",size=19)
  
   df_ggg+ theme_economist(base_size = 28,dkpanel = TRUE,base_family="Avenir")+
     labs(title=expression(paste("Convergence of expected likelihood ",pi)),x="iterations",size=20)
  
plot(df_gg$x, df_gg$y)

Resqfinal_Ut = apply( Resqnew_Ut[,, (maxiter/4) : maxiter],c(1,2), mean )

Resqfinal_U1 <- apply( Resqnew_U1[,, ((maxiter/4):maxiter)],c(1,2),mean)

Resqfinal_Vt <- apply( Resqnew_Vt[,,((maxiter/4) : maxiter)],c(1,2),mean )
Resqfinal_V1 <- apply( Resqnew_V1[,, ((maxiter/4):maxiter)],c(1,2),mean)

b1_ut = mean( mat_b[( (maxiter/4) : maxiter),1] )
bL_ut = mean( mat_b[( (maxiter/4) : maxiter),2] )
b1_u1 = mean( mat_b[( (maxiter/4) : maxiter),3] )
bL_u1 = mean( mat_b[((maxiter/4) : maxiter),4] )
b1_vt = mean( mat_b[((maxiter/4) : maxiter),5] )
bL_vt = mean( mat_b[((maxiter/4) : maxiter),6] )
b1_v1 = mean( mat_b[( (maxiter/4) : maxiter),7] )
bL_v1 = mean( mat_b[((maxiter/4) : maxiter),8] )



