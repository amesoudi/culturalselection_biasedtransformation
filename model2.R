
# load packages------------------------------------
library(emg)
library(RColorBrewer)

# n <- 10  # number of agents per group
# g <- 5  # number of groups
# Beta <- 0.5  # skewness of exponentially modified Gaussian (EMG) / strength of biased transmission (technically the inverse of lambda rate parameter). 0=gaussian, 1=weak, 2=strong
# s <- 0  # probability of selecting the highest value individual as demonstrator, as opposed to picking at random
# Sigma <- 1  # sd of EMG / strength of modification: 0=off, 1=weak, 2=strong
# d <- 0.5  # fitness increment for B, and (doubled) for C lineage, relative to A
# m <- 0.1  # migration or inter-group contact rate
# t_max <- 100  # maximum time steps

philtrans_model2 <- function (n, g, Beta, s, Sigma, d, m, t_max) {
  
  # create a matrix with t_max rows and n*g columns (one per agent), filled with NAs, then convert to data.frame
  output <- as.data.frame(matrix(NA, t_max, n*g))
  names(output) <- paste(paste("g", rep(1:g,each = n), sep = ""), paste("a", rep(1:n,g), sep = ""), sep = "_")  # purely cosmetic: rename the columns with g1_a1, g1_a2 etc. (g=group, a=agent)
  
  Fst <- rep(NA, t_max)  # second output vector for storing Fst, one per generation
  
  # lookup table of traits and their fitnesses
  trait_fitness <- data.frame(trait = c("none","X","A1","A2","A3","A4","A5","B1","B2","B3","B4","B5","C1","C2","C3","C4","C5"), fitness = c(-1,0,1,2,3,4,5,1+d,2+d,3+d,4+d,5+d,1+2*d,2+2*d,3+2*d,4+2*d,5+2*d))
  
  # function to get Fst
  get_Fst <- function(agent) {
    
    g <- max(agent$group)  # no. of subpops; no need for n as all are equal
    traits <- unique(agent$trait)  # list of unique traits in this generation
    s <- length(traits)  # no of unique traits
    
    X <- 0
    for (i in 1:s) {
      X <- X + mean(agent$trait == traits[i])^2  # sum up all the squared means of each variant i
    }
    total.var <- 1 - X  # 1 - sum of squared means
    
    X <- rep(0,g)  # variation per group
    for (j in 1:g) {  # for each group j
      for (i in 1:s) {  # for each variant i
        X[j] <- X[j] + mean(agent$trait[agent$group == j] == traits[i])^2  # sum the squared means of each variant i
      }
      X[j] <- 1 - X[j]  # 1 - sum of squared means within group j
    }
    within.var <- mean(X)  # mean of all variances across all groups
    
    Fst <- (total.var - within.var) / total.var  # return Fst
    
    if (is.na(Fst)) {  # if Fst is NaN due to only 1 variant present, set to 0
      Fst <- 0
    }
    
    Fst  # return Fst
  }
  
  # create agents with group indicator, "none" trait slot and fitness
  agent <- data.frame(group = rep(1:g, each = n), trait = rep("none", n*g), fitness = rep(NA, n*g))  # agent <- data.frame(group = rep(1:g, each = n), trait = sample(trait_fitness$trait, n*g, replace = T), fitness = rep(NA, n*g))  # create first generation
  
  # set fitness based on trait
  agent$fitness <- trait_fitness$fitness[match(agent$trait, trait_fitness$trait)]
  
  output[1,] <- agent$trait  # add first generation's traits to first row
  
  Fst[1] <- get_Fst(agent)  # add Fst for first generation
  
  for (t in 2:t_max) {
    
    # CULTURAL SELECTION
    # copy agent to previous_agent dataframe
    previous_agent <- agent
    
    # vector of whether agents select the highest-value demonstrator (TRUE), otherwise copy at random (FALSE)
    select_best <- runif(n*g) < s
    
    # for each group, agents copy trait of the highest-value or random demonstrator within that group
    for (i in 1:g) {
      
      # for select_best in group i, copy the trait of the agent with the highest fitness in that group (if there are ties, pick randomly)
      agent$trait[agent$group == i & select_best] <- sample(previous_agent$trait[previous_agent$fitness == max(previous_agent$fitness[previous_agent$group == i]) & previous_agent$group == i], 1, replace = T)
      
      # for copy-at-random in group i, pick randomly from traits 
      agent$trait[agent$group == i & !select_best] <- sample(previous_agent$trait[previous_agent$group == i], sum(agent$group == i & !select_best), replace = TRUE)
      
    }
    
    # BIASED TRANSFORMATION / RANDOM MUTATION
    
    change <- remg(n*g, mu = 0, sigma = Sigma, lambda = 1/Beta)  # draw continuous modification
    change[change >= 1] <- 1  # convert to discrete change (up/down the trait lineage, or no change)
    change[change <= -1] <- -1
    change[change > -1 & change < 1] <- 0
    
    for (i in 1:(g*n)) {
      
      # if change is positive/biased, move towards X
      if (change[i] == 1) {
        trait <- agent$trait[i]
        
        if (trait == "none") agent$trait[i] <- "X"  # none goes to X
        if (trait == "A5") agent$trait[i] <- "A4"  # otherwise go down the levels, towards X
        if (trait == "A4") agent$trait[i] <- "A3"
        if (trait == "A3") agent$trait[i] <- "A2"
        if (trait == "A2") agent$trait[i] <- "A1"
        if (trait == "A1") agent$trait[i] <- "X"
        if (trait == "B5") agent$trait[i] <- "B4"
        if (trait == "B4") agent$trait[i] <- "B3"
        if (trait == "B3") agent$trait[i] <- "B2"
        if (trait == "B2") agent$trait[i] <- "B1"
        if (trait == "B1") agent$trait[i] <- "X"
        if (trait == "C5") agent$trait[i] <- "C4"
        if (trait == "C4") agent$trait[i] <- "C3"
        if (trait == "C3") agent$trait[i] <- "C2"
        if (trait == "C2") agent$trait[i] <- "C1"
        if (trait == "C1") agent$trait[i] <- "X"
      }
      
      # if change is negative/against the bias, move away from X
      if (change[i] == -1) {
        
        trait <- agent$trait[i]
        
        if (trait == "X") agent$trait[i] <- sample(c("A1","B1","C1"), 1, prob = c(1/3,1/3,1/3))  # X goes to A1, B1 or C1 with equal probability
        if (trait == "A1") agent$trait[i] <- "A2"   # then go up the levels, away from X
        if (trait == "A2") agent$trait[i] <- "A3"
        if (trait == "A3") agent$trait[i] <- "A4"
        if (trait == "A4") agent$trait[i] <- "A5"
        if (trait == "B1") agent$trait[i] <- "B2"
        if (trait == "B2") agent$trait[i] <- "B3"
        if (trait == "B3") agent$trait[i] <- "B4"
        if (trait == "B4") agent$trait[i] <- "B5"
        if (trait == "C1") agent$trait[i] <- "C2"
        if (trait == "C2") agent$trait[i] <- "C3"
        if (trait == "C3") agent$trait[i] <- "C4"
        if (trait == "C4") agent$trait[i] <- "C5"
      }
      
    }
    
    # MIGRATION
    
    # get list of migrants, based on probability m
    migrants <- runif(n*g) < m  
    # randomly put migrants back into migrant slots, ignoring group (i.e. Wright's island model)
    agent$trait[migrants] <- sample(agent$trait[migrants], length(agent$trait[migrants]), replace = FALSE)
    
    # UPDATE FITNESSES based on new trait
    agent$fitness <- trait_fitness$fitness[match(agent$trait, trait_fitness$trait)]
    
    # STORE RESULTS
    output[t,] <- agent$trait
    Fst[t] <- get_Fst(agent)
    
  }
  
  list(output = output, Fst = Fst)  # export from function
}
  
# run model here---------------------------------------------

m2 <- philtrans_model2(n = 10, g = 3, Beta = 0, s = 1, Sigma = 1, d = 0.5, m = 0, t_max = 100)
View(m2$output)
plot(m2$Fst, type = 'l', ylim = c(0,1), ylab = "Fst", xlab = "timestep")

# representative model dynamics------------------------------

# Beta = 10 (v. strong bias) plus s=0 (no selection) causes convergence on X (actually doesn't need to be so strong)
m2 <- philtrans_model2(n = 100, g = 5, Beta = 10, s = 0, Sigma = 1, d = 0.5, m = 0, t_max = 500)
View(m2$output)
plot(m2$Fst, type = 'l', ylim = c(0,1), ylab = "Fst", xlab = "timestep")

# Beta = 0 (no bias) plus s=1 (strong selection) causes convergence on different A5/B5/C5 values
m2 <- philtrans_model2(n = 100, g = 5, Beta = 0, s = 1, Sigma = 1, d = 0.5, m = 0, t_max = 500)
View(m2$output)
plot(m2$Fst, type = 'l', ylim = c(0,1), ylab = "Fst", xlab = "timestep")

# Beta = 0, s = 1, m = 0.1 (high migration) causes convergence on same trait, not always highest (premature convergence)
m2 <- philtrans_model2(n = 100, g = 5, Beta = 0, s = 1, Sigma = 1, d = 0.5, m = 0.1, t_max = 500)
View(m2$output)
plot(m2$Fst, type = 'l', ylim = c(0,1), ylab = "Fst", xlab = "timestep")

# Beta = 0, s = 1, m = 0.01 (low migration) causes convergence on highest trait (C5)
m2 <- philtrans_model2(n = 100, g = 5, Beta = 0, s = 1, Sigma = 1, d = 0.5, m = 0.01, t_max = 500)
View(m2$output)
plot(m2$Fst, type = 'l', ylim = c(0,1), ylab = "Fst", xlab = "timestep")


# visualise change over time---------------------------

visual <- function (n, g, Beta, s, Sigma, d, m, t_max) {
  
  m2 <- philtrans_model2(n, g, Beta, s, Sigma, d, m, t_max)
  m2 <- m2$output
  m2[m2 == "none"] <- ""  # change "none" to blank
  
  m2_color <- m2  # set colors for each trait type
  m2_color[m2 == ""] <- "#ffffff"
  m2_color[m2 == "X"] <- brewer.pal(n = 9, name = "Greys")[4]
  m2_color[m2 == "A1"] <- brewer.pal(n = 6, name = "Oranges")[1]
  m2_color[m2 == "A2"] <- brewer.pal(n = 6, name = "Oranges")[2]
  m2_color[m2 == "A3"] <- brewer.pal(n = 6, name = "Oranges")[3]
  m2_color[m2 == "A4"] <- brewer.pal(n = 6, name = "Oranges")[4]
  m2_color[m2 == "A5"] <- brewer.pal(n = 6, name = "Oranges")[5]
  m2_color[m2 == "B1"] <- brewer.pal(n = 6, name = "Greens")[1]
  m2_color[m2 == "B2"] <- brewer.pal(n = 6, name = "Greens")[2]
  m2_color[m2 == "B3"] <- brewer.pal(n = 6, name = "Greens")[3]
  m2_color[m2 == "B4"] <- brewer.pal(n = 6, name = "Greens")[4]
  m2_color[m2 == "B5"] <- brewer.pal(n = 6, name = "Greens")[5]
  m2_color[m2 == "C1"] <- brewer.pal(n = 6, name = "Purples")[1]
  m2_color[m2 == "C2"] <- brewer.pal(n = 6, name = "Purples")[2]
  m2_color[m2 == "C3"] <- brewer.pal(n = 6, name = "Purples")[3]
  m2_color[m2 == "C4"] <- brewer.pal(n = 6, name = "Purples")[4]
  m2_color[m2 == "C5"] <- brewer.pal(n = 6, name = "Purples")[5]
  
  par(mgp=c(2,1,0))
  
  plot(NULL, xlim = c(1,nrow(m2)), ylim = c(1,ncol(m2)), ylab = "individual", xlab = "timestep", main = paste("β = ",Beta,", s = ",s, ", m = ",m, sep = ""), yaxt = "n")
  
  for (i in 1:(g-1)) {
    abline(h = ncol(m2)/g*i + 0.5, lty = 2)
  }
  
  for (i in 1:ncol(m2)) { # i denotes agent
    for (j in 1:nrow(m2)) {  # j denotes timestep
      points(j,i, pch = 22, bg = m2_color[j,i], cex = 1.5, lwd = 0.5)
      text(j,i,m2[j,i], cex = 0.4, font = 2)
    }
  }
}

# create visualisation here-----------------------------
visual(n = 5, g = 3, Beta = 2, s = 0.9, Sigma = 1, d = 0.5, m = 0, t_max = 200)

# create fig 4 -------------------------------------

# NB may not be the same every time...

png("fig4.png", width = 30, height = 20, units = 'cm', res = 600)

par(mfrow=c(2,2)) # 2 rows, 2 columns

visual(n = 5, g = 3, Beta = 2, s = 0, Sigma = 1, d = 0.5, m = 0, t_max = 50) # top left, convergence on X via BT
visual(n = 5, g = 3, Beta = 0, s = 0.5, Sigma = 1, d = 0.5, m = 0, t_max = 50) # top right, divergence via CS
visual(n = 5, g = 3, Beta = 0, s = 0.5, Sigma = 1, d = 0.5, m = 0.1, t_max = 50) # bottom left, convergence on C5 via CS+migration

dev.off()

