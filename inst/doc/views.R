## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ffp)
library(ggplot2)
set.seed(321)

# invariance / stationarity
x <- diff(log(EuStockMarkets))

dim(x)

head(x)

## -----------------------------------------------------------------------------
# Returns 20% higher than average
mu <- mean(x[ , "FTSE"]) * 1.2

# ffp views constructor
views <- view_on_mean(x = as.matrix(x[ , "FTSE"]), mean = mu)
views

## -----------------------------------------------------------------------------
# Prior probabilities
p_j <- rep(1 / nrow(x), nrow(x)) 

## -----------------------------------------------------------------------------
# Solve Minimum Relative Entropy (MRE)
ep <- entropy_pooling(
  p      = p_j, 
  Aeq    = views$Aeq, 
  beq    = views$beq, 
  solver = "nlminb"
)
ep

## -----------------------------------------------------------------------------
mu <- c(
  mean(x[ , "FTSE"]) * 1.2, # Return 20% higher than average
  mean(x[ , "CAC"]) * 0.9   # Return 10% bellow average
)

# ffp views constructor
views <- view_on_mean(x = as.matrix(x[ , c("FTSE", "CAC")]), mean = mu)
views

## -----------------------------------------------------------------------------
# Solve MRE
ep <- entropy_pooling(
  p      = p_j, 
  Aeq    = views$Aeq, 
  beq    = views$beq, 
  solver = "nlminb"
)
ep

## -----------------------------------------------------------------------------
class(ep)

## ---- fig.width=7, fig.align='center'-----------------------------------------
autoplot(ep) + 
  scale_color_viridis_c(option = "C", end = 0.75) + 
  labs(title    = "Posterior Probability Distribution", 
       subtitle = "Views on Expected Returns", 
       x        = NULL, 
       y        = NULL)

## -----------------------------------------------------------------------------
conditional   <- ffp_moments(x, ep)$mu 
unconditional <- colMeans(x)

conditional / unconditional - 1

## -----------------------------------------------------------------------------
# opinion
vol_cond <- sd(x[ , "CAC"]) * 0.9

# views
views <- view_on_volatility(x = as.matrix(x[ , "CAC"]), vol = vol_cond)
views

## -----------------------------------------------------------------------------
ep <- entropy_pooling(
  p      = p_j, 
  Aeq    = views$Aeq, 
  beq    = views$beq, 
  solver = "nlminb"
)
ep

## ---- fig.width=7, fig.align='center'-----------------------------------------
autoplot(ep) + 
  scale_color_viridis_c(option = "C", end = 0.75) + 
  labs(title    = "Posterior Probability Distribution", 
       subtitle = "View on Volatility", 
       x        = NULL, 
       y        = NULL)

## -----------------------------------------------------------------------------
unconditional <- apply(x, 2, sd)
conditional   <- sqrt(diag(ffp_moments(x, ep)$sigma))

conditional / unconditional - 1

## -----------------------------------------------------------------------------
cor_view <- cor(x) # unconditional correlation matrix
cor_view["FTSE", "DAX"] <- 0.83 # perturbation (investor belief)
cor_view["DAX", "FTSE"] <- 0.83 # perturbation (investor belief)

## -----------------------------------------------------------------------------
views <- view_on_correlation(x = x, cor = cor_view)
views

## -----------------------------------------------------------------------------
ep <- entropy_pooling(
  p      = p_j, 
  Aeq    = views$Aeq, 
  beq    = views$beq, 
  solver = "nlminb"
)
ep

## ---- fig.width=7, fig.align='center'-----------------------------------------
autoplot(ep) + 
  scale_color_viridis_c(option = "C", end = 0.75) + 
  labs(title    = "Posterior Probability Distribution", 
       subtitle = "View on Correlation", 
       x        = NULL, 
       y        = NULL)

## -----------------------------------------------------------------------------
cov2cor(ffp_moments(x, p = ep)$sigma)

## -----------------------------------------------------------------------------
views <- view_on_rank(x = x, rank = c(2, 1))
views

## -----------------------------------------------------------------------------
ep <- entropy_pooling(
  p      = p_j, 
  A      = views$A, 
  b      = views$b, 
  solver = "nloptr"
)
ep

## ---- fig.width=7, fig.align='center'-----------------------------------------
autoplot(ep) + 
  scale_color_viridis_c(option = "C", end = 0.75) + 
  labs(title    = "Posterior Probability Distribution", 
       subtitle = "View on Ranking/Momentum", 
       x        = NULL, 
       y        = NULL)

## -----------------------------------------------------------------------------
conditional   <- ffp_moments(x, ep)$mu  
unconditional <- colMeans(x)

conditional / unconditional - 1

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(ghyp)

# multivariate t distribution
dist_fit <- fit.tmv(data = x, symmetric = TRUE, silent = TRUE)
dist_fit

## -----------------------------------------------------------------------------
set.seed(123)
# random numbers from `dist_fit`
simul  <- rghyp(n = 100000, object = dist_fit)

## -----------------------------------------------------------------------------
p_z <- rep(1 / 100000, 100000)
views <- view_on_marginal_distribution(
  x     = x, 
  simul = simul, 
  p     = p_z
)
views

## -----------------------------------------------------------------------------
ep <- entropy_pooling(
  p      = p_j, 
  Aeq    = views$Aeq, 
  beq    = views$beq, 
  solver = "nloptr"
)
ep

## ---- fig.width=7, fig.align='center'-----------------------------------------
autoplot(ep) + 
  scale_color_viridis_c(option = "C", end = 0.75) + 
  labs(title    = "Posterior Probability Distribution", 
       subtitle = "View on Marginal Distribution", 
       x        = NULL, 
       y        = NULL)

## -----------------------------------------------------------------------------
ens(ep)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(copula)

# copula (pseudo-observation)
u <- pobs(x)

# copula estimation
cop <- fitCopula(
  copula = claytonCopula(dim = ncol(x)), 
  data   = u, 
  method = "ml"
)

# simulation
r_cop <- rCopula(
  n      = 100000, 
  copula = claytonCopula(param = cop@estimate, dim = 4)
)

## -----------------------------------------------------------------------------
views <- view_on_copula(x = u, simul = r_cop, p = p_z)
views

## -----------------------------------------------------------------------------
ep <- entropy_pooling(
  p      = p_j, 
  Aeq    = views$Aeq, 
  beq    = views$beq, 
  solver = "nloptr"
)
ep

## ---- fig.width=7, fig.align='center'-----------------------------------------
autoplot(ep) + 
  scale_color_viridis_c(option = "C", end = 0.75) + 
  labs(title    = "Posterior Probability Distribution", 
       subtitle = "View on Copulas", 
       x        = NULL, 
       y        = NULL) 

## -----------------------------------------------------------------------------
cop@estimate <- 5

## -----------------------------------------------------------------------------
# simulation
r_cop <- rCopula(
  n      = 100000, 
  copula = claytonCopula(param = cop@estimate, dim = 4)
)

# #views
views <- view_on_copula(x = u, simul = r_cop, p = p_z)

# relative entropy
ep <- entropy_pooling(
  p      = p_j, 
  Aeq    = views$Aeq, 
  beq    = views$beq, 
  solver = "nloptr"
)

## ---- fig.width=7, fig.align='center'-----------------------------------------
autoplot(ep) + 
  scale_color_viridis_c(option = "C", end = 0.75) + 
  labs(title    = "Posterior Probability Distribution", 
       subtitle = "View on Copulas", 
       x        = NULL, 
       y        = NULL) 

