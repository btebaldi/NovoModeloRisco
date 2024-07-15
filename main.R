rm(list =ls())
library(readxl)
library(dplyr)
library(ggplot2)

source(file = "hanssen.R")

# User defined Functions --------------------------------------------------

d.ln <- function(x){
  log(x) - dplyr::lag(log(x))
}

# Data load ---------------------------------------------------------------

IBOV <- read_excel("C:/Users/bteba/Downloads/New folder/IBOV e DOLAR.xlsx", 
                   sheet = "IBOV")
colnames(IBOV) <- c("Data", "Ibov.Abe", "Ibov.Fec")

Dolar <- read_excel("C:/Users/bteba/Downloads/New folder/IBOV e DOLAR.xlsx", 
                    sheet = "DOLAR")
colnames(Dolar) <- c("Data", "Dol.Abe", "Dol.Fec")

tbl <- inner_join(IBOV, Dolar, by = "Data")

rm(list = c("IBOV", "Dolar"))


# Calculo de retornos -----------------------------------------------------

tbl <- tbl %>% 
  mutate(r_bov = d.ln(Ibov.Fec),
         r_dol = d.ln(Dol.Fec)) %>% 
  na.omit()

summary(tbl)

# Parametrizacao t-hanssen (Incondicional) --------------------------------


mopt.Inc <- mle.thanssen(data = tbl$r_bov, lambda = 0, eta = 2.1)

mopt.Inc$par

tbl$r_bov_f1 = dthanssen(x = tbl$r_bov,
                         lambda = mopt.Inc$par["lambda"],
                         eta = mopt.Inc$par["eta"])

ggplot(tbl) + 
  geom_line(aes(x = r_bov, y = r_bov_f1, colour = "t-hanssen"), size = 1) +
  geom_histogram(aes(x = r_bov, y =after_stat(density)), bins = 50, alpha = 0.3, colour = "black") +
  labs()
summary(tbl)

# Parametrizacao t-hanssen (condicional) ----------------------------------

condSample <- (nrow(tbl)-250):nrow(tbl)

mopt.Cond <- mle.thanssen(data = tbl$r_bov[condSample], lambda = 0, eta = 2.1)

mopt.Cond$par

tbl$r_bov_f2 = dthanssen(x = tbl$r_bov, 
                         lambda = mopt.Cond$par["lambda"], 
                         eta = mopt.Cond$par["eta"])

ggplot(tbl) + 
  geom_line(aes(x = r_bov, y = r_bov_f1, colour = "Incond"), size = 1) +
  geom_line(aes(x = r_bov, y = r_bov_f2, colour = "Cond"), size = 1) +
  geom_histogram(aes(x = r_bov, y =after_stat(density)), bins = 50, alpha = 0.3, colour = "black") +
  labs()

ggplot(tbl[condSample,]) + 
  geom_line(aes(x = r_bov, y = r_bov_f1, colour = "Incond"), size = 1) +
  geom_line(aes(x = r_bov, y = r_bov_f2, colour = "Cond"), size = 1) +
  geom_histogram(aes(x = r_bov, y =after_stat(density)), bins = 50, alpha = 0.3, colour = "black") +
  labs()


qthanssen(p=0.9996,
          lambda = mopt.Inc$par["lambda"],
          eta = mopt.Inc$par["eta"], lower.tail = FALSE)


qthanssen(p=0.9996,
          lambda = mopt.Cond$par["lambda"],
          eta = mopt.Cond$par["eta"], lower.tail = FALSE)



# Sorteio de cenarios quantitativos ---------------------------------------

x = rthanssen(n = 10000,
              lambda = mopt.Cond$par["lambda"],
              eta = mopt.Cond$par["eta"] )


MCC <- tibble(x=x,
              f = dthanssen(x = x,
                            lambda = mopt.Cond$par["lambda"],
                            eta = mopt.Cond$par["eta"]))

# KO = 0.25 escolhido a priori
MCC <- MCC %>% filter(abs(x) <= 0.25)

ggplot(MCC) + 
  geom_line(aes(x = x, y = f, colour = "t-hanssen"), size = 1) +
  geom_histogram(aes(x = x, y =after_stat(density)), bins = 100, alpha = 0.3, colour = "black") +
  labs()



# Determinacao de cenarios ------------------------------------------------

cen_hist <- list()

for(i in seq_len(nrow(tbl)-10)){
  M <- list(tbl[i:(i+9),c("r_bov", "r_dol")] %>% data.matrix() %>% t())
  names(M) <- sprintf("CEN_%d", i)
  colnames(M[[1]]) <- paste("HP", 1:10, sep="")
  cen_hist <- append(cen_hist, M)
}


tbl %>% dplyr::select(x = r_bov) %>% bind_rows(MCC) %>% 
  mutate(f = dthanssen(x = x,
                       lambda = mopt.Cond$par["lambda"],
                       eta = mopt.Cond$par["eta"])) %>% 
  ggplot() + 
  geom_line(aes(x = x, y = f, colour = "t-hanssen"), size = 1) +
  geom_histogram(aes(x = x, y =after_stat(density)), bins = 100, alpha = 0.3, colour = "black") +
  labs()




tudo <- tbl %>% dplyr::select(x = r_bov) %>% bind_rows(MCC) %>% 
  mutate(f = dthanssen(x = x,
                       lambda = mopt.Cond$par["lambda"],
                       eta = mopt.Cond$par["eta"]))

quantile(exp(tudo$x)-1, probs = c(0.99, 0.995, 0.9996, (1-0.99), (1-0.995), (1-0.9996))) * sqrt(2)

tail(sort(exp(tudo$x)-1)) *1.4142



