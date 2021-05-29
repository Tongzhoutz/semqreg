
main_model <- function(a,b,maxiter){

  ptm <- proc.time()
for ( iter in 1:maxiter){
  print(iter)

  for ( j in 2:draws){
    
  Matdraws <- matrix(NA,N,T)
  Matdraws[,1] <- Nu_chain1[, j-1]
  Matdraws[,2] <- Nu_chain2[, j-1]
  Matdraws[,3] <- Nu_chain3[, j-1]
  Matdraws[,4] <- Nu_chain4[, j-1]
  Matdraws[,5] <- Nu_chain5[, j-1]
  Matdraws[,6] <- Nu_chain6[, j-1]
  
  
  ## U_1
  
  Matdraws[,1] <- Nu_chain1[,j-1] + sqrt(var_prop1) * rnorm(N,0,1)
  newObj <- posterior_cpp2(Matdraws,DT,Herm_age,regressor_V1,regressor_U1, Mat_init_Vt ,Mat_init_V1 ,Mat_init_Ut,Mat_init_U1,b1_ut, bL_ut, b1_u1,bL_u1, b1_vt,bL_vt,b1_v1,bL_v1 ,a,b,tauList,meanY,stdY )
  r = pmin( ones_vec, newObj/Obj_chain[,j-1]  )
  prob <- runif(N,0,1)
  Obj_chain[,j] <- if_else( prob <= r , newObj, Obj_chain[, j-1]   )
  Nu_chain1[,j] <- if_else(prob <= r, Matdraws[,1], Nu_chain1[, j-1] )
  Matdraws[,1] <- Nu_chain1[,j]
  acc1[,j] <- prob <=r
  ## U_t
  Matdraws[,2] = Nu_chain2[,j-1] + sqrt(var_prop2) * rnorm(N)
  newObj <- posterior_cpp2(Matdraws,DT,Herm_age,regressor_V1,regressor_U1, Mat_init_Vt ,Mat_init_V1 ,Mat_init_Ut,Mat_init_U1,b1_ut, bL_ut, b1_u1,bL_u1, b1_vt,bL_vt,b1_v1,bL_v1 ,a,b,tauList,meanY,stdY )
  r = pmin( ones_vec, newObj/Obj_chain[,j]  )
  prob <- runif(N,0,1)
  Obj_chain[,j] <- if_else( prob <= r , newObj, Obj_chain[, j]   )
  Nu_chain2[,j] <- if_else(prob <= r, Matdraws[,2], Nu_chain2[, j-1] )
  Matdraws[,2] <- Nu_chain2[,j]
  acc2[,j] <- prob <=r
  
  Matdraws[,3] = Nu_chain3[,j-1] + sqrt(var_prop3) * rnorm(N)
  newObj <- posterior_cpp2(Matdraws,DT,Herm_age,regressor_V1,regressor_U1, Mat_init_Vt ,Mat_init_V1 ,Mat_init_Ut,Mat_init_U1,b1_ut, bL_ut, b1_u1,bL_u1, b1_vt,bL_vt,b1_v1,bL_v1 ,a,b,tauList,meanY,stdY )
  r = pmin( ones_vec, newObj/Obj_chain[,j]  )
  prob <- runif(N,0,1)
  Obj_chain[,j] <- if_else( prob <= r , newObj, Obj_chain[, j]   )
  Nu_chain3[,j] <- if_else(prob <= r, Matdraws[,3], Nu_chain3[, j-1] )
  Matdraws[,3] <- Nu_chain3[,j]
  acc3[,j] <- prob <=r
  
  Matdraws[,4] = Nu_chain3[,j-1] + sqrt(var_prop4) * rnorm(N)
  newObj <- posterior_cpp2(Matdraws,DT,Herm_age,regressor_V1,regressor_U1, Mat_init_Vt ,Mat_init_V1 ,Mat_init_Ut,Mat_init_U1,b1_ut, bL_ut, b1_u1,bL_u1, b1_vt,bL_vt,b1_v1,bL_v1 ,a,b,tauList,meanY,stdY )
  r = pmin( ones_vec, newObj/Obj_chain[,j]  )
  prob <- runif(N,0,1)
  Obj_chain[,j] <- if_else( prob <= r , newObj, Obj_chain[, j]   )
  Nu_chain4[,j] <- if_else(prob <= r, Matdraws[,4], Nu_chain4[, j-1] )
  Matdraws[,4] <- Nu_chain4[,j]
  acc4[,j] <- prob <=r
  
  Matdraws[,5] = Nu_chain3[,j-1] + sqrt(var_prop5) * rnorm(N)
  newObj <- posterior_cpp2(Matdraws,DT,Herm_age,regressor_V1,regressor_U1, Mat_init_Vt ,Mat_init_V1 ,Mat_init_Ut,Mat_init_U1,b1_ut, bL_ut, b1_u1,bL_u1, b1_vt,bL_vt,b1_v1,bL_v1 ,a,b,tauList,meanY,stdY )
  r = pmin( ones_vec, newObj/Obj_chain[,j]  )
  prob <- runif(N,0,1)
  Obj_chain[,j] <- if_else( prob <= r , newObj, Obj_chain[, j]   )
  Nu_chain5[,j] <- if_else(prob <= r, Matdraws[,5], Nu_chain5[, j-1] )
  Matdraws[,5] <- Nu_chain5[,j]
  acc5[,j] <- prob <=r
  
  Matdraws[,6] = Nu_chain3[,j-1] + sqrt(var_prop6) * rnorm(N)
  newObj <- posterior_cpp2(Matdraws,DT,Herm_age,regressor_V1,regressor_U1, Mat_init_Vt ,Mat_init_V1 ,Mat_init_Ut,Mat_init_U1,b1_ut, bL_ut, b1_u1,bL_u1, b1_vt,bL_vt,b1_v1,bL_v1 ,a,b,tauList,meanY,stdY )
  r = pmin( ones_vec, newObj/Obj_chain[,j]  )
  prob <- runif(N,0,1)
  Obj_chain[,j] <- if_else( prob <= r , newObj, Obj_chain[, j]   )
  Nu_chain6[,j] <- if_else(prob <= r, Matdraws[,6], Nu_chain6[, j-1] )
  Matdraws[,6] <- Nu_chain6[,j]
  acc6[,j] <- prob <=r
  
  acceptrate1[j] <- mean(acc1[,j])
  acceptrate2[j] <- mean(acc2[,j])
  acceptrate3[j] <- mean(acc3[,j])
  acceptrate4[j] <- mean(acc4[,j])
  acceptrate5[j] <- mean(acc5[,j])
  acceptrate6[j] <- mean(acc6[,j])
  print(j / draws)
  
  }
  
 
  acceptrate1 <- replace_na(acceptrate1,0)
  mean(acceptrate1)
  acceptrate2 <- replace_na(acceptrate2,0)
  mean(acceptrate2)
  acceptrate3 <- replace_na(acceptrate3,0)
  mean(acceptrate3)
  acceptrate4 <- replace_na(acceptrate4,0)
  mean(acceptrate4)
  acceptrate5 <- replace_na(acceptrate5,0)
  mean(acceptrate5)
  acceptrate6 <- replace_na(acceptrate6,0)
  mean(acceptrate6)
  
  last_ut <- cbind( Nu_chain1[,draws],Nu_chain2[,draws],Nu_chain3[,draws],Nu_chain4[,draws],
                    Nu_chain5[,draws],Nu_chain6[,draws])
  
  ##### M-step
  ut_2T <- as.vector( last_ut[,2:T] )
  ut_lag <- as.vector(last_ut[,1:5])
  
  Herm_ut_regressor <- sapply(  (ut_lag - meanY)/stdY , hermPolys_3 ) %>% t
  colnames(Herm_ut_regressor) <- str_c("V",11:14)
  
  DT_ut_M <- data.table(Herm_ut_regressor,Herm_age_2T)
  ut_reg_Mstep <-  DT_ut_M[, .(
                    C1 = V11  * V21,
                    C2 = V11 * V22,
                    C3 = V11 * V23,
                    C4 = V12  * V21,
                    C5 = V12 * V22,
                    C6 = V12 * V23,
                    C7 = V13  * V21,
                    C8 = V13 * V22,
                    C9 = V13 * V23,
                    C10 = V14  * V21,
                    C11 = V14 * V22,
                    C12 = V14 * V23)] %>% as.matrix
                    
    fit_ut <- rq( ut_2T ~ ut_reg_Mstep-1, tau = tauList, method = "fn"  )
    Resqnew_Ut[,,iter] <- coef(fit_ut)
    
    fit_u1 <- rq( last_ut[,1] ~ regressor_U1-1, tau = tauList, method = "fn" )
    Resqnew_U1[,,iter] <- coef(fit_u1)  
    
    DT[, ut_stack := as.vector(last_ut)]
    vt_2T <- DT[, ut_stack := as.vector(last_ut) ][year != 98, .(vt_2T = utoty_nor - ut_stack)] %>% 
             as.matrix
    
    DT_vt_15 <- DT[year != 108, .( vt_15 = utoty_nor - ut_stack  )] %>% as.matrix
    DT_vt_15_reg <- sapply( DT_vt_15, hermPolys_2  ) %>% t
    DT_vt_age <- DT[year != 98, .( Herm.age_1, Herm.age_2 )] %>% as.matrix
    vt_reg_Mstep <- data.table( 
                   C1= DT_vt_15_reg[,1] * DT_vt_age[,1],
                   C2= DT_vt_15_reg[,1] * DT_vt_age[,2],
                   C3= DT_vt_15_reg[,2] * DT_vt_age[,1],
                   C4= DT_vt_15_reg[,2] * DT_vt_age[,2],
                   C5= DT_vt_15_reg[,3] * DT_vt_age[,1],
                   C6= DT_vt_15_reg[,3] * DT_vt_age[,2]) %>% as.matrix
      
    fit_vt <- rq(vt_2T ~ vt_reg_Mstep-1, tau = tauList, method = "fn")          
    Resqnew_Vt[,,iter] <- coef(fit_vt)       
          
    v1_dep <- DT[year == 98,.(utoty_nor - ut_stack)]  %>% as.matrix          
    
    fit_v1 <- rq(v1_dep ~ regressor_V1-1, tau = tauList, method = "fn")
    Resqnew_V1[,,iter] <- coef(fit_v1)
    
    #### Normalization
    
    ##back_ut = mean(Resqnew_Ut[4,,iter])
   ## back_vt = mean(Resqnew_Vt[3,,iter])
    Resqnew_Ut[,,iter] <- Resqnew_Ut[,,iter] - apply( Resqnew_Ut[,,iter],1,mean )  ## demean by 3 rows
    Resqnew_Ut[4,,iter] <- Resqnew_Ut[4,,iter] + 2.0 
    Resqnew_Vt[,,iter] <- Resqnew_Vt[,,iter] - apply( Resqnew_Vt[,,iter],1,mean )  ## demean by 3 rows
    Resqnew_Vt[3,,iter] <- Resqnew_Vt[3,,iter] + 2.0 
    
    Resqnew_Ut[1,,iter] <- Resqnew_Ut[1,,iter] - ( (1 - tauList[11])/bL_ut - tauList[1]/b1_ut )
 
    Resqnew_Vt[1,,iter] <- Resqnew_Vt[1,,iter] - ( (1 - tauList[11])/bL_vt - tauList[1]/b1_vt )
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
         mat_b[iter,8] <- bL_vt
         
         ## mean log-likelihood
         
         mat_lik[iter] <-  mean(log(posterior_cpp2(last_ut,DT,Herm_age,regressor_V1,regressor_U1, Mat_init_Vt ,Mat_init_V1 ,Mat_init_Ut,Mat_init_U1,b1_ut, bL_ut, b1_u1,bL_u1, b1_vt,bL_vt,b1_v1,bL_v1 ,a,b,tauList,meanY,stdY )))
         if (is.na(mat_lik[iter] == TRUE )){
           stop('NA happens')}
         
 Obj_chain <- cbind( Obj_chain[,draws], matrix(0,N, draws-1) )
 Nu_chain1 <- cbind( Nu_chain1[,draws], matrix(0,N, draws -1) )
 Nu_chain2 <- cbind( Nu_chain2[,draws], matrix(0,N, draws -1) )
 Nu_chain3 <- cbind( Nu_chain3[,draws], matrix(0,N, draws -1) )
 Nu_chain4 <- cbind( Nu_chain4[,draws], matrix(0,N, draws -1) )
 Nu_chain5 <- cbind( Nu_chain5[,draws], matrix(0,N, draws -1) )
 Nu_chain6 <- cbind( Nu_chain6[,draws], matrix(0,N, draws -1) )
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
 
         ## 
}  

  return(mat_lik)
  proc.time()-ptm


  
  
  
}