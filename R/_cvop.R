
cvop <- function(y,X) {
       
       n <- ncol(X)
       c <- Variable(n)
 
       quant_loss <- function(u,tau) { 0.5* abs(u) + (tau-0.5)*u    }
       
    result <-    tauList  %>% setNames(str_c("C",1:11)) %>%  map_dfc(function(tau){
         obj <- mean( quant_loss(y -X %*% c,  t = tau) )
         prob <- Problem(Minimize(obj))
         solve(prob, solver="ECOS")$getValue(c) 
         prob_data <- get_problem_data(prob,solver="ECOS")
         if (packageVersion("CVXR") > "0.99-7") {
           ECOS_dims <- ECOS.dims_to_solver_dict(prob_data$data[["dims"]])
         } else {
           ECOS_dims <- prob_data$data[["dims"]]
         }
        
         
         solver_output <- ECOSolveR::ECOS_csolve(c = prob_data$data[["c"]],
                                                 G = prob_data$data[["G"]],
                                                 h = prob_data$data[["h"]],
                                                 dims = ECOS_dims,
                                                 A = prob_data$data[["A"]],
                                                 b = prob_data$data[["b"]])
         if (packageVersion("CVXR") > "0.99-7") {
           direct_soln <- unpack_results(prob, solver_output, prob_data$chain, prob_data$inverse_data)
         } else {
           direct_soln <- unpack_results(prob, "ECOS", solver_output)
         }
         direct_soln$getValue(c) 
       })  %>% 
         as.matrix
  return(result)
  
  
}