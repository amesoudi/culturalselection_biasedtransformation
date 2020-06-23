library(emg)

# n <- 100  # number of agents
# Beta <- 1  # skewness of exponentially modified Gaussian (EMG) / strength of biased transmission (technically the inverse of lambda rate parameter). 0=gaussian, 1=weak, 2=strong
# s <- 1  # probability of selecting the highest value individual as demonstrator, as opposed to picking at random
# Sigma <- 1  # sd of EMG / strength of modification: 0=off, 1=weak, 2=strong
# initial <- 10  # initial mean trait value at start of simulation
# t_max <- 100  # maximum time steps
# r_max <- 10  # number of independent runs

philtrans_model1 <- function (n, initial = 10, Beta, s, Sigma, t_max, r_max) {
  
  output <- as.data.frame(matrix(NA, t_max, r_max))  # create a matrix with t_max rows and r_max columns, filled with NAs, then convert to data.frame
  
  names(output) <- paste("run", 1:r_max, sep="")  # purely cosmetic: rename the columns with run1, run2 etc.
  
  for (r in 1:r_max) {
    
    agent <- data.frame(trait = round(rnorm(n,initial,1)))  # create first generation
    
    output[1,r] <- mean(agent$trait)  # add first generation's mean trait value to first row of column r
    
    for (t in 2:t_max) {
      
      # CULTURAL SELECTION
      # copy agent to previous_agent dataframe
      previous_agent <- agent
      
      # vector of whether agents select the highest-value demonstrator (TRUE), otherwise copy at random (FALSE)
      select_best <- runif(n) < s
      
      # selective copiers adopt highest value trait
      agent$trait[select_best] <- max(previous_agent$trait)
      # others copy at random
      agent$trait[!select_best] <- sample(previous_agent$trait, sum(!select_best), replace = TRUE)
      
      # BIASED TRANSFORMATION / RANDOM MUTATION
      # randomly draw trait values from an exponentially modified Gaussian (EMG) distribution, with gaussian mean of the demonstrator's trait, gaussian sd of Sigma (random change), and exponential rate lambda of inverse Beta (directional change)
      agent$trait <- remg(n, mu = agent$trait, sigma = Sigma, lambda = 1/Beta)  

      # UPPER/LOWER TRAIT BOUNDS
      # impose upper limit of 0-100 for new trait values
      agent$trait[agent$trait > 100] <- 100
      agent$trait[agent$trait < 0] <- 0
      
      # STORE RESULTS
      # get mean trait value and put it into output slot for this generation t and run r
      output[t,r] <- mean(agent$trait) 
    }
  }
  
  # list of parameter values to export
  parameters <- data.frame(par_name = c("n", "beta","s","sigma"), par_value = c(n,Beta,s,Sigma))
  
  list(output, parameters)  # export data from function
}

# time series plot function-------------------------------

timeseries_plot <- function (m1) {
  
  output <- as.data.frame(m1[1])
  parameters <- as.data.frame(m1[2])

  # plot a thick line for the mean p
  plot(rowMeans(output), type = 'l', ylab = "mean trait value", xlab = "timestep", ylim = c(0,100), lwd = 3, main = paste("n = ", parameters$par_value[parameters$par_name == "n"], ", β = ", parameters$par_value[parameters$par_name == "beta"], ", s = ", parameters$par_value[parameters$par_name == "s"], ", σ = ", parameters$par_value[parameters$par_name == "sigma"], sep = ""))
  
  for (r in 1:ncol(output)) {  
    lines(output[,r], type = 'l')  # add lines for each run, up to r_max
  }
  
}

# run model here---------------------------------------------

m1 <- philtrans_model1(n = 1000, Beta = 0.5, s = 0.1, Sigma = 1, t_max = 100, r_max = 5)
timeseries_plot(m1)

# examples--------------------------------------------

# biased transmission / no selection:
timeseries_plot(philtrans_model1(n = 1000, Beta = 1, s = 0, Sigma = 1, t_max = 200, r_max = 5)) # r-shaped increase
timeseries_plot(philtrans_model1(n = 10000, Beta = 1, s = 0, Sigma = 1, t_max = 200, r_max = 5)) # n has no effect
# no bias, selection
timeseries_plot(philtrans_model1(n = 1000, Beta = 0, s = 0.1, Sigma = 1, t_max = 100, r_max = 5)) # s-shaped curve
timeseries_plot(philtrans_model1(n = 100, Beta = 0, s = 0.1, Sigma = 1, t_max = 100, r_max = 5)) # n changes speed of diffusion
timeseries_plot(philtrans_model1(n = 1000, Beta = 0, s = 0.1, Sigma = 1, t_max = 100, r_max = 5)) # n changes speed of diffusion
timeseries_plot(philtrans_model1(n = 10000, Beta = 0, s = 0.1, Sigma = 1, t_max = 100, r_max = 5)) # n changes speed of diffusion

# both biased transmission and cultural selection generate directional cultural change
# but there are differences:
# n has no effect on biased transformation, but does affect cultural selection
# BT shows r-shaped diffusion, CS shows s-shaped (see Henrich 2004; although see Reader)

# both combined:
timeseries_plot(philtrans_model1(n = 1000, Beta = 0.5, s = 0.1, Sigma = 1, t_max = 100, r_max = 5)) 


# visualisations of distributions for Fig 1 -----------------------------------------
Beta <- 0
plot(density(remg(1000000, mu = 0, sigma = 0, lambda = 1/Beta)))

# unbiased mutation (beta = 0)
png("unbiased_mutation.png", width = 10, height = 10, units = 'cm', res = 600)

plot(density(remg(10000000, mu = 0, sigma = 1, lambda = 1/0)), main = "unbiased mutation\n(β = 0)", ylab = "probability density", xlab = "change to demonstrator's trait value", xaxt = "n", xlim = c(-4,4), yaxt = "n")
axis(side = 1, at = seq(-4, 4, by = 2), labels = seq(-4, 4, by = 2))

dev.off()

# biased transformation (beta = 2)
png("biased_transformation.png", width = 10, height = 10, units = 'cm', res = 600)

plot(density(remg(10000000, mu = 0, sigma = 1, lambda = 1/1.5)), main = "biased transformation\n(β = 2)", ylab = "probability density", xlab = "change to demonstrator's trait value", xaxt = "n", xlim = c(-4,6), yaxt = "n")
axis(side = 1, at = seq(-4, 6, by = 2), labels = seq(-4, 6, by = 2))

dev.off()


# plot effect of variables on biased transformation and cultural selection for fig 2----------------------

# measure is timestep95, the number of timesteps taken to reach mean value of 95 (not 100 given mutation)

png("fig2.png", width = 25, height = 15, units = 'cm', res = 600)

par(mfrow=c(2,4)) # 2 rows, 4 columns
length.out <- 20

# label for CS
plot(NULL, xlim = c(0,10), ylim = c(0,10), xaxt = "n", yaxt = "n", ylab = "", xlab = "", axes = F)
text(5,5,"Cultural\nSelection\n(CS)", cex = 2)

timeseries_plot(philtrans_model1(n = 10000, Beta = 0, s = 0.05, Sigma = 1, t_max = 100, r_max = 10)) # s-shaped curve

# effect of s, strength of cultural selection
s_values <- seq(0,1, length.out = length.out)
timestep95 <- rep(NA, length(s_values))
for (i in 1:length(s_values)) {
  x <- as.data.frame(philtrans_model1(n = 1000, Beta = 0, s = s_values[i], Sigma = 1, t_max = 1000, r_max = 10)[1])
  x <- rowMeans(x)
  if (any(x > 95)) {
    timestep95[i] <- min(which(x > 95))
  }
  else {
    timestep95[i] <- NA
  }
}
plot(s_values,timestep95, type = 'l', ylim = c(min(timestep95, na.rm = T), max(timestep95, na.rm = T)), xlab = "s", ylab = "timesteps to reach trait 95", main = "Effect of s with β = 0")

# effect of n, population size, with s > 0
n_values <- seq(100,2000, length.out = length.out)
timestep95 <- rep(NA, length(n_values))
for (i in 1:length(n_values)) {
  x <- as.data.frame(philtrans_model1(n = n_values[i], Beta = 0, s = 0.5, Sigma = 1, t_max = 1000, r_max = 10)[1])
  x <- rowMeans(x)
  if (any(x > 95)) {
    timestep95[i] <- min(which(x > 95))
  }
  else {
    timestep95[i] <- NA
  }
}
plot(n_values,timestep95, type = 'l', ylim = c(max(timestep95, na.rm = T)*0.6, max(timestep95, na.rm = T)*1.1), xlab = "n", ylab = "timesteps to reach trait 95", main = "Effect of n with s = 0.5, β = 0")  # decreases with n

# label for BT
plot(NULL, xlim = c(0,10), ylim = c(0,10), xaxt = "n", yaxt = "n", ylab = "", xlab = "", axes = F)
text(5,5,"Biased\nTransformation\n(BT)", cex = 2)

timeseries_plot(philtrans_model1(n = 10000, Beta = 2, s = 0, Sigma = 1, t_max = 100, r_max = 10)) # r-shaped curve

# effect of Beta, strength of biased transformation
Beta_values <- seq(0,10, length.out = length.out)
timestep95 <- rep(NA, length(Beta_values))
for (i in 1:length(Beta_values)) {
  x <- as.data.frame(philtrans_model1(n = 1000, Beta = Beta_values[i], s = 0, Sigma = 1, t_max = 1000, r_max = 10)[1])
  x <- rowMeans(x)
  if (any(x > 95)) {
    timestep95[i] <- min(which(x > 95))
  }
  else {
    timestep95[i] <- NA
  }
}
plot(Beta_values,timestep95, type = 'l', ylim = c(0,max(timestep95, na.rm = T)), xlab = "β", ylab = "timesteps to reach trait 95", main = "Effect of β with s = 0")

# effect of n, population size, with Beta = 1
n_values <- seq(100,2000, length.out = length.out)
timestep95 <- rep(NA, length(n_values))
for (i in 1:length(n_values)) {
  x <- as.data.frame(philtrans_model1(n = n_values[i], Beta = 2, s = 0, Sigma = 1, t_max = 1000, r_max = 10)[1])
  x <- rowMeans(x)
  if (any(x > 95)) {
    timestep95[i] <- min(which(x > 95))
  }
  else {
    timestep95[i] <- NA
  }
}
plot(n_values,timestep95, type = 'l', ylim = c(0,max(timestep95, na.rm = T)*2), xlab = "n", ylab = "timesteps to reach trait 95", main = "Effect of n with β = 2, s = 0")  # flat, no effect

dev.off()
