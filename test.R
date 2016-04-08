# test.R
#
# the purpose of this script is to test the functions in global.R
#
# first version: 160331
# this version:  160407
# last change by: Yujiao Li

#source("global.R", local = TRUE)

source("SCM_data_simulation.R", local = TRUE)
source("SCM_estimation.R", local = TRUE)
source("SCM_plot.R", local = TRUE)

# set value for test
timeLength = 20
seed.ts = 66
tau = 0.3
controlNumber = 10
omega = 0.3
seed.control = 1 
invTime = 3
alpha = 1
rou = 0.5
seed.obs.random = 66

seed.sample = 10
sample.percent = 0.5




#==========================================    SCM_data_simulation.R  ===============================================#
# Function_1 
ts.data <- TSdataSim(timeLength = timeLength, seed.ts = seed.ts, tau = tau)


# Function_2
mst.data <- MTSdata.Sim( controlNumber = controlNumber, timeLength = timeLength,
             tau = tau, omega = omega, seed.control = seed.control )


# Function_3
y.hat <- ts.data
y.obs <- Obs.Sim(y.hat = y.hat, tau = tau, invTime =invTime ,alpha = alpha)

# Function_4
y.rand <- ObsRand.Sim(y.obs, rou = rou, tau = tau, 
            controlNumber = controlNumber, timeLength = timeLength, seed.obs.random = 3)

# Function_5
data.scm <- DataSCM(controlNumber = controlNumber, timeLength = timeLength , invTime = invTime,
         tau = tau, omega = omega, seed.control = seed.control, alpha = alpha)

# Function_6
sample.data <- SampleData(data.popu = data.scm , seed.sample = 10, sample.percent = 0.5)


#========================================       SCM_estimation.R     =============================================#

# Function_7
SCM.estimate(data = data.scm, invTime = invTime, Y.name = "Treat.Hat")


# Function_8
integ.data <- DataIntegrate.plot(
            timeLength = timeLength,
            controlNumber = controlNumber,
            invTime = invTime,
            
            alpha = alpha,
            tau = tau,
            omega = omega,
            rou = rou,
            
            seed.control = seed.control,
            seed.obs.random = seed.obs.random,
            seed.sample = seed.sample,
            
            sample.percent = sample.percent
)


# Function_9
raw.data <- integ.data$data
Plot.SCM(raw.data, invTime =invTime, randomOption = TRUE)
