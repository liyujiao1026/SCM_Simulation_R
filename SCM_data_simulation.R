# SCM_data_simulation.R script
#
# the purpose of this script is to simulate the Synthetic Control Method
#
# first version: 160331
# this version:  160411
# last change by: Yujiao Li

#==================================================================================================#

library(ggplot2)
library(tidyr)


#=== Funcion_1: TSdataSim() ===========================================================================================#
# TSdata.Sim() generate one time series data:
# Y_t = tau * (Y_t-1) + et
# @timeLength: the length of time series
# @seed.ts: the random seed

TSdata.Sim <- function(timeLength, seed.ts, tau) {
            n <- timeLength
            set.seed(seed.ts)
            e <- rnorm(n, mean = 0, sd = 1)
            y <- c()
            y[1] <- e[1]
            
            for (i in 2:n) {
                        y[i] <- tau * y[i - 1] + e[i]
            }
            
            return(y)
            
}



#=== Function_2: MTSdata.Sim() =====================================================================================================#
# MTSdata.Sim() generate multiple time series data as the potential control units

# Model: Z_i = omega * Z_0 + ( 1- omega) * Z_ii,
#        Z_0 is one common fixed time series inclued in each Zi,
#        Z_ii is the individual's time series
# return: Z_i and one added column of their mean value as the y.hat

# @controlNumber: the number of control units
# @timeLength: the length of time period

MTSdata.Sim <- function(controlNumber, timeLength, tau, omega, seed.control) {
            
            set.seed(seed.control)
            seed.control.series <- sample(1:100000 , controlNumber )
            
            Z0 <- TSdata.Sim(seed = seed.control, tau = tau, timeLength = timeLength)
            Yt <- sapply(seed.control.series, TSdata.Sim, timeLength = timeLength, tau = tau)
            
            Zt <- apply(Yt, 2, function(x) {
                        omega * Z0 + (1 - omega) * x
            })
            
            data <- data.frame(Zt)
            
            rownames(data) <- paste("Time", 1:timeLength , sep = "_")
            colnames(data) <- paste("ControlUnit", 1:controlNumber ,sep = "_")
            
            return(data)
}



#=== Function_3: Obs.Sim() ==================================================================# 

# Obs.Sim  generate the observation 
# Aefore intervention, Treat.Obs = Treat.Hat 
# After intervention,  Treat.Obs = Treat.Hat + alpha
# where alpha is the caused effect from intervention

# return: y.obs

# @ y.hat: the equally weighted average of other control units
# @ intervention_Time : which time point is the intervention time
# @ alpha : the caused effect
# @ seed.post : random walk for the post-intervention observations

Obs.Sim <- function(y.hat, tau, invTime, alpha) {
            
            y.obs <- c()
            # (1) before the intervention, y.hat==y.observation
            y.obs[1:invTime] <- y.hat[1:invTime] 
            
            
            # (2) after the intervention, y.hat = y.obs + alpha
            l <- length(y.hat)
            y.obs[invTime:l] <- y.hat[invTime:l] + alpha
            return(y.obs)
}






#=== Function_4: ObsRand.Sim() =====================================================================# 

# ObsRand.Sim() generate the observation which is random variable
# y.obs.random <- y.obs + rou * epsilo
# E (y.obs.random) = y.obs
# Var(y.obs.random) = var(sum(w^2 * sigma(z)  +  rou^2 * sigma(epsilo)))
# sigma(z) = sqrt(1/(1-tau^2))


# return: y.obs.rand

# @ y.obs: the y observation which is the epectation of each y_t
# @ rou: the weight on the epsilo
# @ seed.obs.random: the random seed for simulation of random y.observation


ObsRand.Sim <- function(y.obs, rou, tau,
                        controlNumber, timeLength, seed.obs.random) {
            
            set.seed(seed.obs.random)
            y.obs.random <- c()
            
            if (tau < 1 & tau > 0) {
                        var_y.obs <- (1 / controlNumber) * (1 / (1 - tau ^ 2))  +  (rou ^ 2)
                        
                        for (i in 1:timeLength) {
                                    y.obs.random[i] <- rnorm(1, mean = y.obs[i], sd = sqrt(var_y.obs))
                        }
                        
            } else if (tau == 1) {     
                        for (k in 1:timeLength) {  
                                    var_y.obs.k <- (1 / controlNumber) * (k)  +  (rou ^ 2)
                                    y.obs.random[k] <- rnorm(1, mean = y.obs[k], sd = sqrt(var_y.obs.k))
                        }
                        
            } else {stop("tau cannot be larger than 1 for stationary time series.")}
            
            return(y.obs.random)
}



#=== Funcion_5: DataSCM() ===============================================================================================#
#
# DataSCM(): generate the dataframe with control units and treated unit
#             
#               will calculate the difference: alpha = y.Hat - y.Obs
#               (2.1) y.Hat ==> Synthetic value (weighted average) and
#               (2.2) y.Obs ==> y.Hat + alpha (after intervention)
#
# return:  dataframe for SCM estimation


DataSCM <- function(timeLength , invTime , controlNumber,
                    tau, omega , alpha, seed.control) {
            
            # 1. simulate the data
            data <- MTSdata.Sim(controlNumber = controlNumber, timeLength = timeLength , 
                                tau = tau, omega = omega, seed.control = seed.control)
            
            y.Hat <- apply(data, 1, mean)
            
            y.Obs <- Obs.Sim(y.hat = y.Hat, invTime = invTime, 
                             alpha = alpha, tau = tau )
            y.Obs.rand <- ObsRand.Sim(y.obs = y.Obs, rou = rou, tau = tau,
                                      controlNumber = controlNumber, timeLength = timeLength, 
                                      seed.obs.random = seed.obs.random )
            
            # 2.Merge as the dataframe
            data.Scm <- data.frame("Time" = 1:timeLength, data, 
                                   "Treat.Hat" = y.Hat, "Treat.Obs" = y.Obs, "Treat.Obs.Rand" = y.Obs.rand) 
            
            return(data.Scm)
            
}



#=== Funcion_6: SampleData() ===============================================================================================#
#
# SampleData(): extract the subset of data generated from DataSCM()
# return:  sample of data 

SampleData <- function(data.popu , seed.sample, sample.percent) {
            set.seed(seed.sample)
            
            control.col <- grep(pattern = "Control", x = colnames(data.popu) , value = FALSE)
            sample.size <- (length(control.col) )* sample.percent
            sample.r0 <- sample(control.col, size = sample.size)
            
            basic.col <- which( !(1:ncol(data.scm))  %in% control.col )
            data <- data.popu[ , sort(c(basic.col, sample.r0)) ]
            
            return(data)
}


