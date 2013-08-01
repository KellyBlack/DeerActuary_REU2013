# ###############################################################
# R script to approximate the balance in an insurance
# fund with a focus on deer collisions.
#
# Authors: Lindong Zhou
#          Kelly Black
#

TESTING         = TRUE   # Will print graphs if true
INITIALIZE_SEED = FALSE  # Will initialize the seed if true

if (TESTING) {
  # The script is set up for testing purposes.
  
  if (INITIALIZE_SEED) {
    # Set the seed using the value in the data file.
    load('randomSeed.Rdata')
    set.seed(theSeed)
  }
  
  # Reset the graphics context
  if (length(dev.list()) > 0) {
    dev.off()
  }
}

T <- 10       # Set the final time
N <- 25000    # Set the number of sample points
r <- 1.5			# r tilde - The scaled growth factor
F <- 25000	  # F tilde - The scaled carrying capacity
alpha <- 0.05  # Noise coeff. in the deer equation.

# Calculate the time step and generate the random numbers
dt <- T/N
t <- seq(0, T, dt)
dw <- sqrt(dt) * rnorm(N)
w <- c(0, cumsum(dw))

# Set the initial values and calculate the scaled parameters.
x0 <- F
z0 <- F/x0
a <- alpha ^ 2 - r
b <- alpha
g0 <- 1 - r/(0.5 * b ^ 2 - a)
phi <- exp((a - 0.5 * b ^ 2) * seq(0, T, dt) - b * w)

# Approximate the stochastic integral for the z equation.
# TODO - convert to a Milstein approximation.
integral <- vector(length = N + 1)
integral <- c(0, cumsum(exp((0.5 * b ^ 2 - a) * t[1:N] + b * w[1:N]) * dw + 
                          0.5 * b * exp((0.5 * b ^ 2 - a) * t[1:N] + 
                                          b * w[1:N]) * (dw * dw - dt)))

# Calculate the values of Z based on the approximation of the integral
# above. Then undo the transformation to calculate the population.
Ztrue <- r/(0.5 * b ^ 2 - a) + g0 * phi - 
  b * r/(0.5 * b ^ 2 - a) * phi * integral
Xtrue <- F/Ztrue

if (TESTING) {
  # Plot the deer population.
  plot(t, Xtrue,
       type = "l", col = "orange", lwd = 2,       
       xlab = "T", ylab = "X", main = "Deer Population")
}

# Approximate the value of the transformed variable used to
# approximate the population.
L <- 5000
R <- N/L
Dt <- R * dt
Xmil <- vector(length = L + 1)
Xmil[1] <- x0
Xtemp <- x0
for (i in 1:L) {
  Winc <- sum(dw[(R * (i - 1) + 1):(R * i)])
  Xtemp <- Xtemp + Dt * r * Xtemp * (1 - Xtemp/F) + alpha * Xtemp * Winc + 
    0.5 * alpha * alpha * Xtemp * (Winc * Winc - Dt)
  Xmil[i + 1] <- Xtemp
}

if (TESTING) {
  # Add these new approximations to the plot.
  courseTime = seq(0, T, Dt);
  points(courseTime, Xmil, col = "green", pch = 4)
}

# Approximate the amount of money in the bond fund
# First set the parameters used for the bond fund equation
rho <- 0.004    # Government bond rate
beta <- 17     # Cost of claim proportional to population
P <- 430000       # Premium
g <- 0.04       # Profit margin

# Now approximate the value of the bond fund.
m0 <- (beta * F - P)/(rho - g)
Mmil <- vector(length = N + 1)
Mmil[1] <- m0
Mtemp <- m0
for (i in 1:N) {
  Mtemp <- Mtemp + dt * (rho * Mtemp - beta * Xtrue[i] + P) -
    dw[i] * beta * Xtrue[i] - 0.5 * beta * alpha * Xtrue[i] * 
    (dw[i] * dw[i] - dt)
  if (Mtemp < 0) {
    Mtemp <- 0
  }
  Mmil[i + 1] <- Mtemp
}

if (TESTING) {
  # Plot the value of the funds
  plot(t, Mmil, type = "p", pch = 4, col = "purple",
       xlab = "T", ylab = "M", main = "Insurance Payout")
}

# Simulate population and payout using different alpha, gamma and P

T <- 10       # Set the final time
N <- 25000    # Set the number of sample points
r <- 1.5  		# r tilde - The scaled growth factor
F <- 25000	  # F tilde - The scaled carrying capacity
dt <- T/N     # Calculate the time step
t <- seq(0, T, dt)

# Set the initial values
x0 <- F
z0 <- F/x0

rho <- 0.004    # Government bond rate
beta <- 17     # Cost of claim proportional to population
g <- 0.04       # Profit margin

approx <- function(alpha, P) {   # Noise coefficient and premium
  
  # Generate the random numbers
  dw <- sqrt(dt) * rnorm(N)
  w <- c(0, cumsum(dw))
  
  # Calculate the scaled parameters.
  a <- alpha ^ 2 - r
  b <- alpha
  g0 <- 1 - r/(0.5 * b ^ 2 - a)
  phi <- exp((a - 0.5 * b ^ 2) * seq(0, T, dt) - b * w)
  
  # Approximate the stochastic integral for the z equation.
  # TODO - convert to a Milstein approximation.
  integral <- vector(length = N + 1)
  integral <- c(0, cumsum(exp((0.5 * b ^ 2 - a) * t[1:N] + b * w[1:N]) * dw +
                            0.5 * b * exp((0.5 * b ^ 2 - a) * t[1:N] + 
                                            b * w[1:N]) * (dw * dw - dt)))
  
  # Calculate the values of Z based on the approximation of the integral
  # above. Then undo the transformation to calculate the population.
  Ztrue <- r/(0.5 * b ^ 2 - a) + g0 * phi - 
    b * r/(0.5 * b ^ 2 - a) * phi * integral
  Xtrue <- F/Ztrue
  
  # Now approximate the value of the bond fund.
  m0 <- (beta * F - P)/(rho - g)
  Mmil <- vector(length = N + 1)
  Mmil[1] <- m0
  Mtemp <- m0
  for (i in 1:N) {
    Mtemp <- Mtemp + dt * (rho * Mtemp - beta * Xtrue[i] + P) -
      dw[i] * beta * alpha * Xtrue[i] -
      0.5 * beta * alpha * alpha * Xtrue[i] * (dw[i] * dw[i] - dt)
    Mmil[i + 1] <- Mtemp
  }
  
  # Return the values of X and M at time T
  return (c(Xtrue[N + 1], Mmil[N + 1]))
  
}

# Create the file and write out a header
setwd("C:/Users/Linton/Desktop")
cat(file = "deerSim1.csv", "alpha, P, X, M\n", append=FALSE)

for (i in 1:1000) {
  App <- approx(0.05, 475000)
  cat(file = "deerSim1.csv", 0.05, 475000, App[1], App[2], append = TRUE, 
      sep = ",", fill = TRUE)
}

deerSim1 <- read.csv(file = "deerSim1.csv", header = TRUE)
hist(deerSim1$X, probability = TRUE, xlab = "X", ylab = "Frequency density",
     main = "Distribution of X", col = "grey")
qqnorm(deerSim1$X)
qqline(deerSim1$X)
hist(deerSim1$M, probability = TRUE, xlab = "M", ylab = "Frequency density",
     main = "Distribution of M", col = "grey")
qqnorm(deerSim1$M)
qqline(deerSim1$M)

# Create the file and write out a header
setwd("C:/Users/Linton/Desktop")
cat(file = "deerSim2.csv", "alpha, P, Xmean, Xvar, Mmean, Mvar\n", 
    append=FALSE)

for (i in seq(0.01, 0.1, 0.01)) {
  for (j in seq(430000, 525000, 5000)) {
    Xapp <- Mapp <- vector(length = 1000)
    for (k in 1:1000) {
      App <- approx(i, j)
      Xapp[k] <- App[1]
      Mapp[k] <- App[2]
    }
    cat(file = "deerSim2.csv", i, j, mean(Xapp), var(Xapp), 
        mean(Mapp), var(Mapp), append = TRUE, sep = ",", fill = TRUE)
  }
}