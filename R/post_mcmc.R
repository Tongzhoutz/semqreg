

   post_mcmc   <- function(Mat_draw){
  
      
     ########## 1: density of the Data
      
    
     ## Vectorize 
     Mat_1 <- Mat_draw[,1] %>% as.matrix
     Mat_26 <- Mat_draw[, 2:6] %>% as.matrix
     Mat_15 <- Mat_draw[,1:5] %>% as.matrix
     
     
     
     ##### 1.1: V_t | V_{t-1}
     
     DT_copy <- data.table( DT[, .(year,utoty_nor,age,uc_resid, Herm.age_1,Herm.age_2,Herm.age_3,Herm.age_4)])
     cols <- c("herm_Vect.1", "herm_Vect.2" ,"herm_Vect.3")
     DT_copy[, `:=`(Vect =  utoty_nor - as.vector(Mat_draw) , Mat_draw = as.vector(Mat_draw))  ] 
     Herm_Vect <- DT_copy[, sapply( (Vect-meanY)/stdY, hermPolys_2  )  ] %>% t %>% as.matrix
   
    regressor_Vt_updated <- data.table( C1 = Herm_Vect[1 : 3960,1] * DT[(N+1) : 4752, Herm.age_1],
                                         C2 = Herm_Vect[1 : 3960,1] * DT[(N+1)  : 4752, Herm.age_2],
                                         
                                         C3 = Herm_Vect[1 : 3960,2] * DT[(N+1)  : 4752, Herm.age_1],
                                         C4 = Herm_Vect[1 : 3960,2] * DT[(N+1)  : 4752, Herm.age_2],
                                       
                                         C5 = Herm_Vect[1 : 3960,3] * DT[(N+1)  : 4752, Herm.age_1],
                                        C6 = Herm_Vect[1 : 3960,3] * DT[(N+1)  : 4752, Herm.age_2])  %>% 
                              as.matrix
     
#     regressor_Vt_updated <- data.table( C1 = Herm_Vect[1 : 3960,1] * DT[(N+1) : 4752, Herm.age_1],
#                                         C2 = Herm_Vect[1 : 3960,1] * DT[(N+1)  : 4752, Herm.age_2],
#                                         C4 = Herm_Vect[1 : 3960,2] * DT[(N+1)  : 4752, Herm.age_1],
#                                         C5 = Herm_Vect[1 : 3960,2] * DT[(N+1)  : 4752, Herm.age_2],
#                                         C6 = Herm_Vect[1 : 3960,2] * DT[(N+1)  : 4752, Herm.age_3],
                                         
 #                                        C7 = Herm_Vect[1 : 3960,3] * DT[(N+1)  : 4752, Herm.age_1],
#                                         C8 = Herm_Vect[1 : 3960,3] * DT[(N+1)  : 4752, Herm.age_2],
                                         
#                                         C9 = Herm_Vect[1 : 3960,3] * DT[(N+1)  : 4752, Herm.age_3]
#                                         )  %>% 
#                              as.matrix
     
    predict_vt_updated <- regressor_Vt_updated %*% Mat_init_Vt
    
   DT_copy2 <- data.table( Vect = DT_copy[year != 98, Vect], predict_vt_updated ) %>% 
                set_names( c( "Vect", str_c("tau.",1:11) ) )
   dens2_vt <- DT_copy2[,  `:=`(
      prod12 =  tau_diff/(tau.2 - tau.1) * fifelse( Vect > tau.1 & Vect <= tau.2,1,0 ),
      prod23 =  tau_diff/(tau.3 - tau.2) * fifelse( Vect > tau.2 & Vect <= tau.3,1,0 ),
      prod34 =  tau_diff/(tau.4 - tau.3) * fifelse( Vect > tau.3 & Vect <= tau.4,1,0 ),
      prod45 =  tau_diff/(tau.5 - tau.4) * fifelse( Vect > tau.4 & Vect <= tau.5,1,0 ),
      prod56 =  tau_diff/(tau.6 - tau.5) * fifelse( Vect > tau.5 & Vect <= tau.6,1,0 ),
      prod67 =  tau_diff/(tau.7 - tau.6) * fifelse( Vect > tau.6 & Vect <= tau.7,1,0 ),
      prod78 =  tau_diff/(tau.8 - tau.7) * fifelse( Vect > tau.7 & Vect <= tau.8,1,0 ),
      prod89 =  tau_diff/(tau.9 - tau.8) * fifelse( Vect > tau.8 & Vect <= tau.9,1,0 ),
      prod910 =  tau_diff/(tau.10 - tau.9) * fifelse( Vect > tau.9 & Vect <= tau.10,1,0 ),
      prod1011 =  tau_diff/(tau.11 - tau.10) * fifelse( Vect > tau.10 & Vect <= tau.11,1,0 ),
      prod01 = tauList[1] * b1_vt * exp( b1_vt*(Vect - tau.1) ) * fifelse( Vect <= tau.1,1,0),
      prod11 = (1 - tauList[11]) * bL_vt * exp(-bL_vt * (Vect-tau.11) ) * fifelse(Vect > tau.11,1,0))][, 
                                          .(dens_tt = prod01+prod12+prod23+prod34+prod45+prod56+prod67+prod78+prod89+ prod910 + prod1011 + prod11)][,
                                        .(dens2_vt = dens_tt[1:N]*dens_tt[(N+1):(2*N)]*dens_tt[(2*N+1):(3*N)]*dens_tt[(3*N+1):(4*N)]*dens_tt[(4*N+1):(5*N)])]
      #### 1.2: V_1 | AGE_1                                                                                                                                                                                                                ]
     predict_v1 <- regressor_V1  %*% Mat_init_V1
     DT_copy_v1 <- data.table( Vect = DT_copy[year == 98, Vect], predict_v1 ) %>% 
                   set_names( c("Vect", str_c("tau.",1:11)))
  
     dens2_v1 <-  DT_copy_v1[, `:=`( 
       prod01 = tauList[1] * b1_v1 * exp(b1_v1*(Vect - tau.1)  ) * fifelse(Vect <= tau.1,1,0),
       prod12 =  tau_diff/(tau.2-tau.1) * fifelse( Vect  > tau.1 & Vect  <= tau.2,1,0 ), 
       prod23 =  tau_diff/(tau.3-tau.2) * fifelse( Vect  > tau.2 & Vect  <= tau.3,1,0 ), 
       prod34 =  tau_diff/(tau.4-tau.3) * fifelse( Vect  > tau.3 & Vect  <= tau.4,1,0 ),
       prod45 =  tau_diff/(tau.5-tau.4) * fifelse( Vect  > tau.4 & Vect <= tau.5,1,0 ),
       prod56 =  tau_diff/(tau.6-tau.5) * fifelse( Vect  > tau.5 & Vect  <= tau.6,1,0 ),
       prod67 =  tau_diff/(tau.7-tau.6) * fifelse( Vect  > tau.6 & Vect  <= tau.7,1,0 ),
       prod78 =  tau_diff/(tau.8-tau.7) * fifelse( Vect  > tau.7 & Vect  <= tau.8,1,0 ),
       prod89 =  tau_diff/(tau.9-tau.8) * fifelse(Vect > tau.8 & Vect  <= tau.9,1,0 ),
       prod910 = tau_diff/(tau.10-tau.9) * fifelse( Vect  > tau.9 & Vect  <= tau.10,1,0 ), 
       prod1011 =tau_diff/(tau.11-tau.10) * fifelse( Vect  > tau.10 & Vect  <= tau.11,1,0 ) , 
       prod11 = (1 - tauList[11]) * bL_v1 * exp( -bL_v1*(Vect -tau.11) )*fifelse(Vect  > tau.11,1,0)
     )][,.(dens2_v1 = prod01 + prod12 + prod23 + prod34 + prod45 + prod56 + prod67 + 
             prod78 + prod89 + prod910 + prod1011 + prod11)] 
     
     
     ##### 2: prior, U_t | U_{t-1}, AGE_t
     
     regressor_herm_ut <-  t( sapply( ( as.vector(Mat_15 )-meanY)/stdY, hermPolys_3 ))
     DT_copy_Ut <- DT_copy[year != 98, .(Herm.age_1,Herm.age_2,Herm.age_3)] %>%
       tibble::add_column(regressor_herm_ut) %>% 
        do.call(data.frame,.) %>% 
       as.data.table %>% 
       set_names( c(str_c("Herm.age_",1:3) , str_c("V",1:4)) )
     
     regressor_Ut_updated <- DT_copy_Ut[,
                                   .( 
                                  C1 =V1 * Herm.age_1 ,
                                  C2 =V1 * Herm.age_2,
                                  C3 =V1 * Herm.age_3,      
                                  C4 =V2 * Herm.age_1 ,
                                  C5 =V2 * Herm.age_2,
                                  C6 =V2 * Herm.age_3,
                                  C7 =V3 * Herm.age_1 ,
                                  C8 =V3 * Herm.age_2,
                                  C9 =V3 * Herm.age_3,
                                  C10 =V4 * Herm.age_1 ,
                                  C11 =V4 * Herm.age_2,
                                  C12 =V4 * Herm.age_3
                                  )] %>% 
                             as.matrix
     
     predict_ut_updated <- regressor_Ut_updated %*% Mat_init_Ut
     
     DT_copy_Ut <- data.table( Mat_26 = as.vector(Mat_26) , predict_ut_updated ) %>% 
                  set_names( c("Mat_26", str_c("tau.",1:11)))
     
     dens1_ut <- DT_copy_Ut[,`:=`(
       ifelse01 =  fifelse( Mat_26 <= tau.1, 1,0),
       ifelse12 =  fifelse( Mat_26 > tau.1 & Mat_26 <= tau.2, 1, 0 )  ,
       ifelse23 =  fifelse( Mat_26 > tau.2 & Mat_26 <= tau.3, 1, 0 )  ,
       ifelse34 =  fifelse( Mat_26 > tau.3 & Mat_26 <= tau.4, 1, 0 )  ,
       ifelse45 =  fifelse( Mat_26 > tau.4 & Mat_26 <= tau.5, 1, 0 )  ,
       ifelse56 =  fifelse( Mat_26 > tau.5 & Mat_26 <= tau.6, 1, 0 )  ,
       ifelse67 =  fifelse( Mat_26 > tau.6 & Mat_26 <= tau.7, 1, 0 )  ,
       ifelse78 =  fifelse( Mat_26 > tau.7 & Mat_26 <= tau.8, 1, 0 )  ,
       ifelse89 =  fifelse( Mat_26 > tau.8 & Mat_26 <= tau.9, 1, 0 )  ,
       ifelse910 =  fifelse( Mat_26 > tau.9 & Mat_26 <= tau.10, 1, 0 )  ,
       ifelse1011 =  fifelse( Mat_26 > tau.10 & Mat_26 <= tau.11, 1, 0 ) ,
       ifelse11 =  fifelse(Mat_26 > tau.11, 1, 0),
       tau12 =  tau_diff/(tau.2-tau.1),
       tau23 =  tau_diff/(tau.3-tau.2),
       tau34 =  tau_diff/(tau.4-tau.3),
       tau45 =  tau_diff/(tau.5-tau.4),
       tau56 =  tau_diff/(tau.6-tau.5),
       tau67 =  tau_diff/(tau.7-tau.6),
       tau78 =  tau_diff/(tau.8-tau.7),
       tau89 =  tau_diff/(tau.9-tau.8),
       tau910 =  tau_diff/(tau.10-tau.9),
       tau1011 =  tau_diff/(tau.11-tau.10))][,
                                            `:=`(prod01 = tauList[1] * b1_ut * exp(b1_ut*(Mat_26 - tau.1)  ) * ifelse01,
                                                  prod12 = ifelse12 * tau12, 
                                                  prod23 = ifelse23 * tau23, 
                                                  prod34 = ifelse34 * tau34, 
                                                  prod45 = ifelse45 * tau45, 
                                                  prod56 = ifelse56 * tau56, 
                                                  prod67 = ifelse67 * tau67, 
                                                  prod78 = ifelse78 * tau78, 
                                                  prod89 = ifelse89 * tau89, 
                                                  prod910 = ifelse910 * tau910, 
                                                  prod1011 = ifelse1011 * tau1011, 
                                                  prod11 = (1 - tauList[11]) * bL_ut * exp( -bL_ut*(Mat_26-tau.11) )*ifelse11)][, `:=`( prod01 = if_else(is.na(prod01),0,prod01 ),  prod12 = if_else(is.na(prod12),0,prod12 ), prod23 = if_else(is.na(prod23),0,prod23 ), prod34 = if_else(is.na(prod34),0,prod34 ), prod45 = if_else(is.na(prod45),0,prod45 ), prod56 = if_else(is.na(prod56),0,prod56 ), prod67 = if_else(is.na(prod67),0,prod67 ), prod78 = if_else(is.na(prod78),0,prod78 ),
 prod89 = if_else(is.na(prod89),0,prod89 ), prod910 = if_else(is.na(prod910),0,prod910 ), prod1011 = if_else(is.na(prod1011),0,prod1011 ) ,  prod11 = if_else(is.na(prod11),0,prod11)) ][, 
                                                .(dens_tt = prod01 + prod12 + prod23 + prod34 + prod45 + prod56 + prod67 + prod78 + prod89 + prod910 + prod1011 + prod11)][,.(dens1_ut = dens_tt[1:N]*dens_tt[(N+1):(2*N)]*dens_tt[(2*N+1):(3*N)]*dens_tt[(3*N+1):(4*N)]*dens_tt[(4*N+1):(5*N)])]
   
        ##### 3: prior U_1
     predict_u1 <- regressor_U1 %*% Mat_init_U1
     DT_copy_u1 <- data.table(Mat_1 = as.vector(Mat_1), predict_u1 ) %>% 
                   set_names( c("Mat_1", str_c("tau.",1:11)) )
     
     dens1_u1 <- DT_copy_u1[,`:=`( 
       prod01 = tauList[1] * b1_u1 * exp(b1_u1*(Mat_1 - tau.1)  ) * fifelse(Mat_1 <= tau.1,1,0),
       prod12 =  tau_diff/(tau.2-tau.1) * fifelse( Mat_1 > tau.1 & Mat_1 <= tau.2,1,0 ), 
       prod23 =  tau_diff/(tau.3-tau.2) * fifelse( Mat_1 > tau.2 & Mat_1 <= tau.3,1,0 ), 
       prod34 =  tau_diff/(tau.4-tau.3) * fifelse( Mat_1 > tau.3 & Mat_1 <= tau.4,1,0 ),
       prod45 =  tau_diff/(tau.5-tau.4) * fifelse( Mat_1 > tau.4 & Mat_1 <= tau.5,1,0 ),
       prod56 =  tau_diff/(tau.6-tau.5) * fifelse( Mat_1 > tau.5 & Mat_1 <= tau.6,1,0 ),
       prod67 =  tau_diff/(tau.7-tau.6) * fifelse( Mat_1 > tau.6 & Mat_1 <= tau.7,1,0 ),
       prod78 =  tau_diff/(tau.8-tau.7) * fifelse( Mat_1 > tau.7 & Mat_1 <= tau.8,1,0 ),
       prod89 =  tau_diff/(tau.9-tau.8) * fifelse( Mat_1 > tau.8 & Mat_1 <= tau.9,1,0 ),
       prod910 = tau_diff/(tau.10-tau.9) * fifelse( Mat_1 > tau.9 & Mat_1 <= tau.10,1,0 ), 
       prod1011 =tau_diff/(tau.11-tau.10) * fifelse( Mat_1 > tau.10 & Mat_1 <= tau.11,1,0 ) , 
       prod11 = (1 - tauList[11]) * bL_u1 * exp( -bL_u1*(Mat_1-tau.11) )*fifelse(Mat_1 > tau.11,1,0)
     )][,.(dens1_u1 = prod01 + prod12 + prod23 + prod34 + prod45 + prod56 + prod67 + 
             prod78 + prod89 + prod910 + prod1011 + prod11)] 
    
   
     #### 4: consumption
        dens_c =     DT_copy[, .( dens_con =1/( 1 + exp(uc_resid - a - b*Mat_draw -c*Vect))  )] [, 
                                .( dens_c = dens_con[1:N]*dens_con[ (N+1):(2*N)] * dens_con[(2*N+1):(3*N)]*
                                     dens_con[(3*N+1):(4*N)]*dens_con[(4*N+1):(5*N)]*dens_con[(5*N+1):(6*N)]  ) ]               
      #return( as.matrix(dens1_ut) )                      
     return(   as.matrix(  dens2_v1 * dens2_vt * dens1_u1 * dens1_ut *dens_c  )              )
        
}