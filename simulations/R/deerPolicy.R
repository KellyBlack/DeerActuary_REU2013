# ###############################################################
# R script to approximate the balance in an insurance
# fund with a focus on deer collisions.
#
# Authors: Lindong Zhou
#          Kelly Black
#

TESTING         = TRUE   # Will print graphs if true
INITIALIZE_SEED = FALSE  # Will initialize the seed if true

if(TESTING)
  {
    # The script is set up for testing purposes.

    if(INITIALIZE_SEED)
      {
        # Set the seed using the value in the data file.
        load('randomSeed.Rdata')
        set.seed(theSeed)
      }

    # Reset the graphics context
    if(length(dev.list())>0)
      {
        dev.off()
      }
    par(mfrow=c(2,2))
  }

T <- 10        # Set the final time
N <- 3650      # Set the number of sample points
r <- 10			   # r tilde - The scaled growth factor
F <- 1000000	 # F tilde - The scaled carrying capacity
alpha <- 0.1   # Noise coeff. in the deer equation.

# Calculate the time step and generate the random numbers.
dt <- T/N
dw <- sqrt(dt) * rnorm(N)
w  <- c(0, cumsum(dw))

# Set the initial values and calculate the scaled parameters.
x0  <- F
z0  <- F/x0
a   <- alpha ^ 2 - r
b   <- alpha
g0  <- 1 - r/(0.5 * b ^ 2 - a)
t   <- seq(0,T,dt);
phi <- exp((a - 0.5 * b ^ 2) * t - b * w)

# Approximate the stochastic integral for the z equation.
# TODO - convert to a Milstein approximation.
integral <- vector(length = N + 1)
integral[1] <- 0
integral = c(0,cumsum(exp((0.5 * b ^ 2 - a) * t[1:N] + b * w[1:N]) * dw[1:N]));

# Calculate the values of Z based on the approximation of the integral
# above. Then undo the transformation to calculate the population.
Ztrue <- r/(0.5 * b ^ 2 - a) +
  g0 * phi - b * r/(0.5 * b ^ 2 - a) * phi * integral
Xtrue <- F/Ztrue


if(TESTING)
  {
    # Plot the deer population.
    plot(t, Xtrue,
         type = "l", col = "orange", lwd = 2,
         ylim = c(900000, 1100000),
         xlab = "T", ylab = "X", main = "Deer Population")
  }

# Approximate the value of the transformed variable used to
# approximate the population.
L <- 730
R <- N/L
Dt <- R * dt
Xem <- Xmil <- vector(length = L + 1)
Xem[1] <- Xmil[1] <- x0
Xtemp1 <- Xtemp2 <- x0
for (i in 1:L)
{
  Winc <- sum(dw[(R * (i - 1) + 1):(R * i)])
  Xtemp1 <- Xtemp1 + Dt * r * Xtemp1 * (1 - Xtemp1/F) + alpha * Xtemp1 * Winc
  Xem[i + 1] <- Xtemp1
  Xtemp2 <- Xtemp2 + Dt * r * Xtemp2 * (1 - Xtemp2/F) + alpha * Xtemp2 * Winc + 0.5 * alpha * alpha * Xtemp2 * (Winc * Winc - Dt)
  Xmil[i + 1] <- Xtemp2
}

if(TESTING)
  {
    # Add these new approximations to the plot.
    courseTime = seq(0, T, Dt);
    points(courseTime, Xem, col = "green", pch = 16)
    points(courseTime, Xmil, col = "red", pch = 4, lwd = 2)
    legend("bottomleft", legend = c("Euler", "Milstein"),
           col = c("green", "red"), pch = c(16, 4), lwd = c(1, 2), cex = 0.8)
  }


# Now transform the value of z into the population Distribution
Z10 <- vector(length = 1000)
for (i in 1:1000)
{
  dt <- T/N
  dw <- sqrt(dt) * rnorm(N)
  w <- c(0, cumsum(dw[1:(N - 1)]))
  x0 <- F
  z0 <- F/x0
  a <- alpha ^ 2 - r
  b <- alpha
  g0 <- 1 - r/(0.5 * b ^ 2 - a)
  phi <- exp((a - 0.5 * b ^ 2) * T - b * w[N])
  integral <- 0
  ds <- dt
  s = seq(0, T - ds, ds)
  integral = sum(exp((0.5 * b ^ 2 - a) * s + b * w) * dw)
  Ztrue <- r/(0.5 * b ^ 2 - a) +
    g0 * phi - b * r/(0.5 * b ^ 2 - a) * phi * integral
  Z10[i] <- Ztrue
}

if(TESTING)
  {
    cat("The sample mean of Z10 is ",mean(Z10),"\n",
        "The sample standard deviation of Z10 is ",sd(Z10),"\n");
    hist(Z10, probability = TRUE, xlab = "Z", ylab = "Frequency density",
         main = "Distribution of Z", col = "grey")
    qqnorm(Z10)
  }


# Approximate the amount of money in the bond fund
# First set the parameters used for the bond fund equation
I     <- 0.0005
beta  <- 0.1
P     <- 10000
gamma <- 0.1
g     <- 0.0001

# Now approximate the value of the bond fund.
m0     <- (beta * F - P)/(I - g)
Mmil <- vector(length = N + 1)
Mmil[1] <- m0
Mtemp <- m0
for (i in 1:N)
{
  Mtemp <- Mtemp + dt * (I * Mtemp - beta * Xtrue[i] + P) -
            dw[i] * gamma * Xtrue[i] +
            0.5 * gamma * gamma * Xtrue[i] * (dw[i] * dw[i] - dt)
  Mmil[i + 1] <- Mtemp
}


if(TESTING)
{
  # Plot the value of the funds
  plot(t, Mmil, type = "p", pch = 4, col = "purple",
       xlab = "T", ylab = "M", main = "Insurance Payout")
}