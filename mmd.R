library('np')
library('parallel')
library('pracma')
library('HDCD')
# Simulation Sample
set.seed(63130)
N <- 10**6 # Population size
mc <- 100      # Number of MC iteration
ns <- c(200, 800, 1500)       # Sample size

# Population {{{
## Truth
gamma <- c(1, 1, -1)
# v <- rnorm(n, 0, 2)
v <- runif(N, -2, 3) # try uniform distribution

## Observed
v1 = ceiling(v)
v0 <- v1 - 1
# x <- rnorm(n, 1, 4)
x <- runif(N, 0, 5) # try uniform distribution

## Error
eps <- rnorm(N, 0, 1)

## y
y <- gamma[1] + gamma[2]*v + gamma[3]*x + eps

## pop
pop <- data.frame(y, x, v0, v1)
# }}}

# Monte Carlo Simulation {{{
mc_results <- data.frame()
for (n in ns) {
  for (i in 1:mc) {
    idx <- sample(N, n, replace=FALSE)
    yn <- pop[idx, 'y']
    xn <- pop[idx, 'x']
    v0n <- pop[idx, 'v0']
    v1n <- pop[idx, 'v1']

    X.eval <- data.frame(v0n, v1n, xn)

    # Estimation
    ## estimation of eta -- Nadaraya–Watson estimator
    ### estimate the bandwidth first
    bw <- np::npregbw(yn ~ v0n + v1n + xn, regtype='lc') 

    ### estimate eta Kernel regression
    model <- npreg(bws = bw, gradients = TRUE)  
    eta <- as.matrix(predict(model, newdata=X.eval))

    ## Objective function Q = mean(Q0 + Q1) (see Cerquera et al., 2014 page 4)
    Q <- function (x, v0, v1, params, eta) {
      # params <- gamma
      f1 <- params[1]*1 + params[2]*v1 + params[3]*x
      #Q1 <- sapply(f1, hausdorff_dist, Q=eta) * (f1 < eta)
      # Q1 <- (f1 - eta)**2 * (f1 < eta)
      Q1 <- abs(f1 - eta) * (f1 < eta)
      #Q1 <- hausdorff(f1, eta) * (f1 < eta)
      f0 <- params[1]*1 + params[2]*v0 + params[3]*x
      #Q0 <- sapply(f0, hausdorff_dist, Q=eta) * (f0 > eta)
      # Q0 <- (f0 - eta)**2 * (f0 > eta)
      Q0 <- abs(f0 - eta) * (f0 > eta)
      #Q0 <- hausdorff(f0, eta) * (f0 > eta)
      return(mean(Q1 + Q0))
    }

    ## Optimization
    Q_min <- optim(par=c(0, 0, 0), fn=Q, eta=eta, x=xn, v0=v0n, v1=v1n,
      # lower=c(-2, 0, -2), upper=c(2, 2, 2), 
      method='BFGS',
      control=list(maxit=500, trace=6))

    # Print Estimation results
    print(round(Q_min$par, 2))

    ## Q(true) <= min(Q) + epsilon_N
    epsilon_N <- log(n)/n
    Q_true <- Q(xn, v0n, v1n, gamma, eta)
    (Q_true < (Q_min$value + epsilon_N))
    print(paste0('Q_true: ', round(Q_true, 5), ' Q_min: ', round(Q_min$value, 2)))

    ### Sweep the Reals
    candidates <- expand.grid(
      lapply(Q_min$par, FUN=function(x){seq(x - 1, x + 1, 0.01)})
    )

    # Get set C which Q(C) =< min(Q)
    cl <- makeCluster(detectCores())
    clusterExport(cl, varlist=ls())
    C_star <- candidates[
      parallel::parApply(cl=cl, candidates, 
        FUN=function(c){(Q(xn, v0n, v1n, c, eta) <= (Q_min$value + epsilon_N))}, 
        MARGIN=1)
      , ]
    stopCluster(cl)

    temp <- cbind.data.frame(
      sapply(C_star, max),
      sapply(C_star, min),
      var = c('intercept', 'x', 'v'),
      mc = i,
      n = n
    )
    mc_results(mc_results, temp)
  }
}
# }}}

### Graph -- holding other coef fixed {{{
beta <- expand.grid(
 lapply(Q_min$par, FUN=function(x){seq(x - 2, x + 2, 0.2)})
)
Qt <- numeric(length(beta))
varlist <- c('intercept', 'v', 'x')
par(mfrow=c(1, 3))
for (j in 1:3) {
  gt <- gamma
  for (i in 1:length(beta)) {
    gt[j] <- beta[i]
    Qt[i] <- Q(x, v0, v1, gt, eta)
  }
  plot(x=beta, y=Qt, main=paste0('beta_', varlist[j]), ylim=c(0, 5))
  abline(h=(Q_min$value + epsilon_N), col='red', lty=3)
  abline(v=gamma[j], lty=2)
  abline(v=Q_min$par[j], lty=2, col='blue')
}
### }}}


summary(C_star)




# Another Approach -- Find Vertices

## Population 
get_Cstar <- function(candidates, x, v0, v1, eta) {
  # Function to check if Q1(c, eta) == 0 and Q0(c, eta) == 0
  is_in_Cstar <- function(c) {
    params <- as.numeric(c)
    f1 <- params[1]*1 + params[2]*v1 + params[3]*x
    Q1 <- (f1 - eta)**2 * (f1 < eta)
    
    f0 <- params[1]*1 + params[2]*v0 + params[3]*x
    Q0 <- (f0 - eta)**2 * (f0 > eta)
    
    # Check if both Q1 and Q0 are zero
    all(Q1 == 0) && all(Q0 == 0)
  }
  
  # Apply the condition to all candidates
  Cstar <- candidates[apply(candidates, 1, is_in_Cstar), ]
  return(Cstar)
}

# Example usage
candidates <- expand.grid(
  lapply(Q_min$par, FUN=function(x) { seq(x - 2, x + 2, 0.01) })
)
summary(candidates)

C_star <- get_Cstar(candidates, x, v0, v1, eta)
dim(C_star)
head(C_star)
summary(C_star)

## Sample analogy
get_Cstar <- function(candidates, x, v0, v1, eta, epsilon_N) {
  # Function to check if the constraints are satisfied
  is_in_Cstar <- function(c) {
    params <- as.numeric(c)
    
    # Compute f(x, v0, c) and f(x, v1, c)
    f0 <- params[1]*1 + params[2]*v0 + params[3]*x
    f1 <- params[1]*1 + params[2]*v1 + params[3]*x
    
    # Check the constraints for all n
    all((f0 <= (eta + epsilon_N)) & ((eta - epsilon_N) <= f1))
  }
  
  # Apply the condition to all candidates
  Cstar <- candidates[apply(candidates, 1, is_in_Cstar), ]
  return(Cstar)
}

# Example usage
candidates <- expand.grid(
  lapply(Q_min$par, FUN=function(x) { seq(x - 1, x + 1, 0.01) })
)

# Compute C_star using the new constraints
C_star <- get_Cstar(candidates, x, v0, v1, eta, epsilon_N)
print(C_star)
summary(C_star)

plt1 <- C_star[C_star$Var1 == unique(C_star$Var1)[1], ]

# Plot points with Var1 as x and Var2 as y
par(mfrow=c(1,1))
plot(C_star$Var2, C_star$Var3, 
     xlab = "V", 
     ylab = "X", 
     # main = "Plot of Var1 vs Var2", 
     pch = 19, # Solid circle for points
     col = "blue") # Color of points


plt2 <- C_star[C_star$Var3 == unique(C_star$Var3)[2], ]

# Plot points with Var1 as x and Var2 as y
plot(plt2$Var1, plt2$Var2, 
     xlab = "Var1", 
     ylab = "Var2", 
     main = "Plot of Var1 vs Var2", 
     pch = 19, # Solid circle for points
     col = "blue") # Color of points



