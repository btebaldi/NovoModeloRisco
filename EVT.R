
# Setup -------------------------------------------------------------------
library(readxl)
library(data.table)
library(xts)
library(rugarch)
library(evd)

library(ggplot2)
library(tidyr)
library(dplyr)

rm(list=ls())

# User defined function ---------------------------------------------------

d.ln <- function(x){
  return(log(x) - dplyr::lag(log(x)))
}

determine_KO <- function(cen_hist){
  KO <- rep(NA, 10)
  names(KO) <- sprintf("KO_%d", 1:10)
  for(i in 1:10){
    KO[i] <- quantile(abs(cen_hist[i,1, ]), probs = c(0.9999) )
  }
  return(KO)
}

apply_KO <- function(simulated_data, KO){
  for(i in 10:1){
    removeCol.idx <- which(abs(simulated_data[i, ]) > KO[i])
    simulated_data <- simulated_data[, -removeCol.idx]
  }
  return(simulated_data)
}


# Data load ---------------------------------------------------------------

Ibov <- read_excel("database/DOLAR-IBOV.xlsx", sheet = "IBOV") |> as.data.table()
Dol <- read_excel("database/DOLAR-IBOV.xlsx", sheet = "DOLAR") |> as.data.table()

# Data regularization -----------------------------------------------------

colnames(Ibov) <- c("Date", "Ibov.open", "Ibov.close")
colnames(Dol) <- c("Date", "Dol.open", "Dol.close")

tbl <- merge(Ibov, Dol, by.x = c("Date"), by.y = c("Date"))

# Calcula o log retorno da distribuição
tbl[, c("dl.Ibov", "dl.Dol") := log(.SD) - shift(log(.SD), 1, NA, "lag"), .SDcols=c("Ibov.close", "Dol.close") ]

# remove a primeira linha que contem NA
tbl <- tbl[ -1, .(Date, dl.Ibov, dl.Dol)]



# Construindo cenarios Historicos -----------------------------------------


cen_hist <- array(NA,
                  dim = c(10, 2, nrow(tbl)-10),
                  dimnames = list(sprintf("HP%02d", 1:10),
                                  c("Ibov", "Dol"),
                                  sprintf("Cen_%d", 1:(nrow(tbl)-10))));  

for(i in seq_len(nrow(tbl)-10)){
  cen_hist[ , , i] <- tbl[i:(i+9),c("dl.Ibov", "dl.Dol")] %>% data.matrix() %>% apply(2, cumsum)
  
  # names(M) <- sprintf("CEN_%d", i)
  # colnames(M[[1]]) <- paste("HP", 1:10, sep="")
  # cen_hist <- append(cen_hist, M)
}


# ACF e PACF --------------------------------------------------------------

png(file=sprintf("./figs/%s_ACF_raw.png", "dl.Ibov"), width=6, height=4, units="in", res=100)
TSA::acf(tbl[ , dl.Ibov])
TSA::acf(tbl[ , dl.Ibov^2])
dev.off()

png(file=sprintf("./figs/%s_PACF_raw.png", "dl.Ibov"), width=6, height=4, units="in", res=100)
pacf(tbl[ , dl.Ibov])
dev.off()


# Modelando o Garch -------------------------------------------------------

spec <- ugarchspec( variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(1,0), include.mean = FALSE),
                    distribution.model = "sstd" # Skewed Standardized Student's t-distribution
)

col <- "dl.Ibov"

# times series object
z <- xts(x = tbl[ , get(col)], order.by = as.Date(tbl[ , Date]))
names(z) <- col

# Fit the GARCH model
fit <- ugarchfit(spec = spec, data = z)

residuals <- residuals(fit)/sigma(fit)

# Analise de residuo
png(file=sprintf("./figs/%s_ACF_garch.png", col), width=6, height=4, units="in", res=100)
TSA::acf(residuals)
TSA::acf(residuals^2)
dev.off()

png(file=sprintf("./figs/%s_PACF_garch.png", col), width=6, height=4, units="in", res=100)
pacf(residuals)
dev.off()


# Modelando o EVT ---------------------------------------------------------
hist(residuals, breaks = "FD")


# Install and load necessary packages
# install.packages("skewt")
library(sn)


# Define the log-likelihood function
log_likelihood <- function(params, data) {
  xi <- params[1]
  omega <- params[2]
  alpha <- params[3]
  nu <- params[4]
  
  -sum(sn::dst(data,
               xi = xi,
               omega = omega,
               alpha = alpha, 
               nu = nu,
               log = TRUE))
}

# Initial parameter guesses
initial_params <- c(xi = 0, omega = 1, alpha = 1, nu = 10)

# Perform MLE
mle_fit <- optim(initial_params, log_likelihood, data = as.numeric(residuals),
                 method = "L-BFGS-B",
                 lower = c(-Inf, 0.001, -Inf, 0.001), upper = c(Inf, Inf, Inf, Inf))

# Print the MLE estimates
mle_params <- mle_fit$par
# names(mle_params) <- c("location", "scale", "shape", "df")
names(mle_params) <- c("xi", "omega", "alpha", "nu")
print(mle_params)

# data <- rst(n, location, scale, shape, df)

xx <- sort(as.numeric(residuals))
hist(as.numeric(residuals), breaks = "FD", freq = FALSE)
lines(xx,
      dst(xx,
          xi = mle_params["xi"],
          omega = mle_params["omega"],
          alpha = mle_params["alpha"], 
          nu = mle_params["nu"]), col = "red", lwd = 2)


# Number of samples to generate
n <- 1e4

sim_residuals <- rst(n,
                     xi = mle_params["xi"],
                     omega = mle_params["omega"],
                     alpha = mle_params["alpha"], 
                     nu = mle_params["nu"])



# Use the simulated residuals to generate new time series data Specify the same
# GARCH model used before
spec_sim <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 0), include.mean = FALSE),
  distribution.model = "sstd", # Skewed Standardized Student's t-distribution
  fixed.pars=as.list(coef(fit)) )

# Create an empty matrix to store the simulated data

simulated_garch <- ugarchpath(spec_sim,
                              n.sim = 10,    # The simulation horizon.
                              n.start = 250,  # The burn-in sample.
                              m.sim = 1e4,      # The number of simulations.
                              # presigma = fit@fit$sigma,
                              # prereturns = fit@fit$residuals,
                              preresiduals = sim_residuals)

simulated_data <- fitted(simulated_garch)
colnames(simulated_data) <- sprintf("CEN_%d", 1:ncol(simulated_data))

simulated_data <- apply(simulated_data, MARGIN = 2, cumsum)

png(file=sprintf("./figs/%s_hist_garch_sim.png", col), width=6, height=4, units="in", res=100)
hist(simulated_data[2, ], breaks = "FD")
dev.off()

q_val <- c(0.0004, 0.0050, 0.0100, 0.9900, 0.9950, 0.9996)
dim(simulated_data)
quantile(simulated_data[2, ], probs = q_val)

quantile(cen_hist[2,1, ], probs = q_val)

quantile(tbl$dl.Ibov, probs = q_val)


KO_vector <- determine_KO(cen_hist)

dim(simulated_data)
# simulated_data <- apply_KO(simulated_data, KO = KO_vector)
dim(simulated_data)


png(file=sprintf("./figs/%s_hist_garch_sim with KO.png", col), width=6, height=4, units="in", res=100)
hist(simulated_data[2, ], breaks = "FD")
dev.off()

q_val <- c(0.0004, 0.0050, 0.0100, 0.9900, 0.9950, 0.9996)
dim(simulated_data)
quantile(exp(simulated_data[2, ])-1, probs = q_val)
quantile(exp(cen_hist[2,1, ])-1, probs = q_val)

quantile(tbl$dl.Ibov, probs = q_val)
quantile(cen_hist[1,1, ], probs = q_val)








# Using Conditional information for garch --------------------------------

# Fit the GARCH model
fit.cond <- ugarchfit(spec = spec, data = z["2022/"])

# Analise de residuo

png(file=sprintf("./figs/%s_ACF_Garch_cond.png", col), width=6, height=4, units="in", res=100)
TSA::acf(as.numeric(residuals(fit.cond)))
TSA::acf(as.numeric(residuals(fit.cond))^2)

dev.off()

png(file=sprintf("./figs/%s_PACF_Garch_cond.png", col), width=6, height=4, units="in", res=100)
pacf(as.numeric(residuals(fit.cond)))
dev.off()


# Extract standardized residuals
residuals <- residuals(fit.cond)/sigma(fit.cond)
TSA::acf(residuals)
TSA::acf(residuals^2)

# Modelando o EVT (conditional) -------------------------------------------

hist(residuals, breaks = "FD")


# Initial parameter guesses
initial_params <- c(xi = 0, omega = 1, alpha = 1, nu = 10)

# Perform MLE
mle_fit <- optim(initial_params, log_likelihood, data = as.numeric(residuals),
                 method = "L-BFGS-B",
                 lower = c(-Inf, 0.001, -Inf, 0.001), upper = c(Inf, Inf, Inf, Inf))

# Print the MLE estimates
mle_params <- mle_fit$par
# names(mle_params) <- c("location", "scale", "shape", "df")
names(mle_params) <- c("xi", "omega", "alpha", "nu")
print(mle_params)

# data <- rst(n, location, scale, shape, df)

xx <- sort(as.numeric(residuals))
hist(as.numeric(residuals), breaks = "FD", freq = FALSE)
lines(xx,
      dst(xx,
          xi = mle_params["xi"],
          omega = mle_params["omega"],
          alpha = mle_params["alpha"], 
          nu = mle_params["nu"]), col = "red", lwd = 2)


# Number of samples to generate
n <- 1e4

sim_residuals <- rst(n,
                     xi = mle_params["xi"],
                     omega = mle_params["omega"],
                     alpha = mle_params["alpha"], 
                     nu = mle_params["nu"])



# Use the simulated residuals to generate new time series data Specify the same
# GARCH model used before
spec_sim <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 0), include.mean = FALSE),
  distribution.model = "sstd", # Skewed Standardized Student's t-distribution
  fixed.pars=as.list(coef(fit.cond)) )

# Create an empty matrix to store the simulated data

simulated_garch2 <- ugarchpath(spec_sim,
                               n.sim = 10,    # The simulation horizon.
                               n.start = 250,  # The burn-in sample.
                               m.sim = 1e4,      # The number of simulations.
                               # presigma = fit@fit$sigma,
                               # prereturns = fit@fit$residuals,
                               preresiduals = as.numeric(sim_residuals) )

simulated_data2 <- fitted(simulated_garch2)
colnames(simulated_data2) <- sprintf("CEN_%d", 1:ncol(simulated_data2))


q_val <- c(0.0004, 0.0050, 0.0100, 0.9900, 0.9950, 0.9996)
dim(simulated_data2)
quantile(simulated_data2[2, ], probs = q_val)
quantile(simulated_data[2, ], probs = q_val)

quantile(cen_hist[1,1, ], probs = q_val)
quantile(tbl$dl.Ibov, probs = q_val)

simulated_data11 <- apply_KO(simulated_data, KO = KO_vector)
simulated_data21 <- apply_KO(simulated_data2, KO = KO_vector)
dim(simulated_data11)
dim(simulated_data21)



quantile(simulated_data[2, ], probs = q_val)
quantile(simulated_data2[2, ], probs = q_val)
quantile(cen_hist[2,1, ], probs = q_val)


quantile(simulated_data11[2, ], probs = q_val)
quantile(simulated_data21[2, ], probs = q_val)
quantile(cen_hist[2,1, ], probs = q_val)

quantile(c(simulated_data11[2, ], cen_hist[2,1, ]), probs = q_val)
quantile(c(simulated_data21[2, ], cen_hist[2,1, ]), probs = q_val)






