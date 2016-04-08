# SCM_estimation.R script

# the purpose of this script is to estimate the parameter and 
# generate the data for plot in shiny app
#
# first version: 160331
# this version:  160411
# last change by: Yujiao Li

#==================================================================================================#

library(kernlab)
library(LowRankQP)


#=== Funcion_7: SCM.estimate() ===============================================================================================#

# SCM.estimate(): generate the dataframe with control units and treated unit
# return : According to the dataset and the fitted responsevariable to return the following 4 results:
# (1) w_scm.hat: weight estimated from SCM for each control unit
# (2) Y_scm.hat: based on new weight.hat
# (3) MSE_scm_before: mean square error for the pre-intervention
# (4) alpha1_scm.hat: estimated alpha 

# @ Y.name is the character name of the target variable to estimate ("Treat.Obs" or "Treat.Obs.Rand")

SCM.estimate <- function(data, invTime, Y.name,
                       margin.ipop = 0.0005,
                       sigf.ipop = 5,
                       bound.ipop = 10,
                       quadopt = "ipop") {
          
            
            control.col <- grep(pattern = "Control", colnames(data) , value = FALSE)
            Z <- data[,control.col] %>% as.matrix #control
            
            treat.col <- which(colnames(data) == Y.name)
            Y.obs <- data[,treat.col]  %>% as.matrix    #treat
            
            
            #pre-intervention of Z_matrix and Y_obs were used for estimation 
            Z_matrix <- Z[1:invTime,] %>% as.matrix
            Y_obs <- Y.obs[1:invTime,] %>% as.matrix
            
            
            # set up QP problem
            H <- t(Z_matrix) %*%  (Z_matrix) 
            a <- Y_obs
            c <- -1*c(t(a) %*%  (Z_matrix) )
            A <- t(rep(1, length(c)))
            b <- 1
            l <- rep(0, length(c))
            u <- rep(1, length(c))
            r <- 0
            
            # run QP and obtain w weights
            # ipop
            if (quadopt == "ipop") {
                        
                        res <- ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r, bound = bound.ipop,
                                    margin = margin.ipop, maxiter = 1000, sigf = sigf.ipop)
                        
                        solution.w <- as.matrix(primal(res))
                        
                        
            } else {
                        
                        # LowRankQP
                        if (quadopt == "LowRankQP") {
                                    
                                    res <- LowRankQP(Vmat = H,dvec = c,Amat = A,
                                                     bvec = 1,uvec = rep(1,length(c)),method = "LU")
                                    solution.w <- as.matrix(res$alpha)
                        } 
            }
            
            
            
            # results    
            
            MSE_pre.inv <- as.numeric(t(Y_obs - Z_matrix %*% solution.w) %*% (Y_obs - Z_matrix %*% solution.w )) / nrow(Z_matrix)
            Y.predict <- Z %*% solution.w
            
            # calculate alpha 1 time after intervention
            t.1 <- invTime + 1
            Y_predict.1 <- (as.vector(Z %*% solution.w))[t.1] 
            Y_observe.1 <- as.vector(data$Treat.Obs)[t.1] 
            alphaHat.1 <- (Y_observe.1 - Y_predict.1)
            
            
            return(list("w_scm.hat" = solution.w, "MSE_scm_before" = MSE_pre.inv, 
                        "Y_scm.hat" = Y.predict, "alpha1_scm.hat" = alphaHat.1 ))
}




#=== Funcion_8: DataIntegrate() ===============================================================================================#

# DataIntegrate() is to estimate the new alphaHat seperately assumed of fix observation and random observation

# return: 
# (1)Treat.Hat  
# (2)Treat.Obs (3)Treat.Hat_SCM 
# (4)Treat.Obs.random (5) Treat.Hat_SCM.Rand


DataIntegrate.plot <- function(timeLength = timeLength,
                               controlNumber = controlNumber,
                               invTime = invTime,
                               
                               alpha = alpha,
                               tau = tau,
                               omega = omega,
                               rou = rou,
                               
                               seed.control = seed.control,
                               seed.obs.random = seed.obs.random,
                               seed.sample = seed.sample,
                               
                               sample.percent = sample.percent){
            
            
            # (1) generate the data
            data.scm <- DataSCM(controlNumber = controlNumber, timeLength = timeLength , invTime = invTime,
                                tau = tau, omega = omega, seed.control = seed.control, alpha = alpha)
            
            sample.data <- SampleData(data.popu = data.scm , seed.sample = seed.sample, sample.percent = sample.percent)
            
            # (2).data with fixed & random observation to estimate alpha and Y.hat
            # sampling the data and estimated the alphaHat when observation are fixed
            
            fixObs.result <- SCM.estimate(data = sample.data, invTime = invTime, Y.name = "Treat.Obs")
            randObs.result <- SCM.estimate(data = sample.data, invTime = invTime, Y.name = "Treat.Obs.Rand")
            
            
            sample.data$Treat.Hat_SCM <- fixObs.result$Y_scm.hat
            sample.data$Treat.Hat_SCM.Rand <- randObs.result$Y_scm.hat

            
            # output the dataset and alpha estimated
            data <- data.frame(sample.data)
            
            alphaHat.fix <- fixObs.result$alpha1_scm.hat
            alphaHat.rand <- randObs.result$alpha1_scm.hat
            
            
            return(list( "data" = data, "alpha.fix" = alphaHat.fix, "alpha.rand" = alphaHat.rand))
            }
























