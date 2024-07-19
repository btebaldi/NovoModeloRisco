
# Setup -------------------------------------------------------------------
library(readxl)
library(data.table)
library(xts)
library(rugarch)
library(evd)


# User defined function ---------------------------------------------------

d.ln <- function(x){
  return(log(x) - dplyr::lag(log(x)))
}

determine_KO <- function(simulated_data){
  KO <- rep(NA, 10)
  names(KO) <- sprintf("KO_%d", 1:10)
  for(i in seq_along(KO)){
    KO[i] <- quantile(abs(simulated_data[i, ]), probs = c(0.995) )
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


# ACF e PACF --------------------------------------------------------------

TSA::acf(tbl[ , dl.Ibov])
TSA::acf(tbl[ , dl.Dol])

pacf(tbl[ , dl.Ibov])
pacf(tbl[ , dl.Dol])

# Modelando o Garch -------------------------------------------------------

spec <- ugarchspec( variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(10,0), include.mean = FALSE),
                    distribution.model = "sstd" # Skewed Standardized Student's t-distribution
)

col <- "dl.Ibov"

# times series object
z <- xts(x = tbl[ , get(col)], order.by = as.Date(tbl[ , Date]))
names(z) <- col

# Fit the GARCH model
fit <- ugarchfit(spec = spec, data = z)

# Analise de residuo
TSA::acf(as.numeric(residuals(fit)))
pacf(as.numeric(residuals(fit)))

# z <- xts(x = sp500ret, order.by = as.Date(rownames(sp500ret)))
# garchspec <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
#                         variance.model = list(model ="gjrGARCH"),
#                         distribution.model = "std")
# 
# garchroll <- ugarchroll(garchspec, data = z, n.start = 2500,
#                         refit.window = "moving", refit.every = 100)


# Extract standardized residuals
residuals <- residuals(fit, standardize = TRUE)

# Modelando o EVT ---------------------------------------------------------

# Fit a generalized extreme value (GEV) distribution to the residuals
gev_fit <- evd::fgev(residuals)

# Print the fitted parameters
print(gev_fit)


# Generate random samples from the fitted GEV distribution

# Number of samples to generate
n <- 1e6

sim_residuals <- evd::rgev(n, 
                           loc = gev_fit$estimate["loc"],
                           scale = gev_fit$estimate["scale"],
                           shape = gev_fit$estimate["shape"])

hist(sim_residuals, breaks = "FD")

# Use the simulated residuals to generate new time series data Specify the same
# GARCH model used before
spec_sim <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(10, 0), include.mean = FALSE),
  distribution.model = "sstd", # Skewed Standardized Student's t-distribution
  fixed.pars=as.list(coef(fit)) )

# Create an empty matrix to store the simulated data

simulated_garch <- ugarchpath(spec_sim,
                              n.sim = 10,    # The simulation horizon.
                              n.start = 250,  # The burn-in sample.
                              m.sim = 100000,      # The number of simulations.
                              # presigma = fit@fit$sigma,
                              # prereturns = fit@fit$residuals,
                              preresiduals = sim_residuals,
                              fixed.pars=as.list(coef(fit)))

simulated_data <- fitted(simulated_garch)
colnames(simulated_data) <- sprintf("CEN_%d", 1:ncol(simulated_data))

simulated_data <- apply(simulated_data, MARGIN = 2, cumsum)

hist(simulated_data[2, ], breaks = "FD")

KO_vector <- determine_KO(simulated_data)

dim(simulated_data)
simulated_data <- apply_KO(simulated_data, KO = KO_vector)
dim(simulated_data)
hist(simulated_data[2, ], breaks = "FD")

# Using Conditional information for garch --------------------------------

# Fit the GARCH model
fit.cond <- ugarchfit(spec = spec, data = z["2021/"])

# Analise de residuo
TSA::acf(as.numeric(residuals(fit.cond)))
pacf(as.numeric(residuals(fit.cond)))

# Extract standardized residuals
residuals.cond <- residuals(fit.cond, standardize = TRUE)

# Modelando o EVT (conditional) -------------------------------------------

# Fit a generalized extreme value (GEV) distribution to the residuals
gev_fit.cond <- evd::fgev(residuals.cond)

# Print the fitted parameters
print(gev_fit.cond)

# Generate random samples from the fitted GEV distribution

# Number of samples to generate
n <- 1e6

sim_residuals.cond <- evd::rgev(n, 
                           loc = gev_fit.cond$estimate["loc"],
                           scale = gev_fit.cond$estimate["scale"],
                           shape = gev_fit.cond$estimate["shape"])

hist(sim_residuals.cond, breaks = "FD")

# Use the simulated residuals to generate new time series data Specify the same
# GARCH model used before
spec_sim.cond <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(10, 0), include.mean = FALSE),
  distribution.model = "sstd", # Skewed Standardized Student's t-distribution
  fixed.pars=as.list(coef(fit.cond)) )

# Create an empty matrix to store the simulated data
simulated_garch.cond <- ugarchpath(spec_sim.cond,
                              n.sim = 10,    # The simulation horizon.
                              n.start = 250,  # The burn-in sample.
                              m.sim = 100000,      # The number of simulations.
                              # presigma = fit@fit$sigma,
                              # prereturns = fit@fit$residuals,
                              preresiduals = sim_residuals,
                              fixed.pars=as.list(coef(fit.cond)))

simulated_data.cond <- fitted(simulated_garch.cond)
colnames(simulated_data.cond) <- sprintf("CEN_%d", 1:ncol(simulated_data.cond))

simulated_data.cond <- apply(simulated_data.cond, MARGIN = 2, cumsum)

hist(simulated_data.cond[2, ], breaks = "FD")

dim(simulated_data.cond)
KO_vector <- determine_KO(simulated_data.cond)
simulated_data.cond <- apply_KO(simulated_data.cond, KO = KO_vector)
dim(simulated_data.cond)

hist(simulated_data.cond[2, ], breaks = "FD")


# Graficos ----------------------------------------------------------------



simulated_data.cond %>% 
  as_tibble() %>% 
  mutate(HP = 1:10) %>% 
  pivot_longer(cols = -HP) %>% 
  ggplot() + 
  geom_line(aes(x = HP, y = value, colour = name)) +
  theme_bw() +
  theme(legend.position = "none") + 
  labs()


simulated_data %>% 
  as_tibble() %>% 
  mutate(HP = 1:10) %>% 
  pivot_longer(cols = -HP) %>% 
  ggplot() + 
  geom_line(aes(x = HP, y = value, colour = name)) +
  theme_bw() +
  theme(legend.position = "none") + 
  labs()


ggplot() + 
  geom_density(aes(x = value, fill = "Historical"), colour = "black", alpha=0.5, data = as_tibble(simulated_data[ 2, ])) +
  geom_density(aes(x = value, fill = "Conditional"), colour = "black", alpha=0.5, data = as_tibble(simulated_data.cond[ 2, ])) +
  theme_bw() +
  theme(legend.position = "bottom") + 
  labs(fill = NULL)

q = c(0.99, 0.995, 0.996)
q = sort(c(1-q, q))

quantile(simulated_data[2,], prob = q)
quantile(simulated_data.cond[2,], prob = q)
