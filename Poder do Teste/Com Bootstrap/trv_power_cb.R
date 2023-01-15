
library(munsell)
library(tidyverse)
library(magrittr)
library(reshape2)
library(progress)
knitr::opts_chunk$set(echo = FALSE)

newton.raphson <- function(f, a, b, tol = 1e-5, n = 1000) {
  require(numDeriv) # Package for computing f'(x)
  
  x0 <- a # Set start value to supplied lower bound
  k <- n # Initialize for iteration results
  
  # Check the upper and lower bounds to see if approximations result in 0
  fa <- f(a)
  if (fa == 0.0) {
    return(a)
  }
  
  fb <- f(b)
  if (fb == 0.0) {
    return(b)
  }
  
  for (i in 1:n) {
    dx <- genD(func = f, x = x0)$D[1] # First-order derivative f'(x0)
    x1 <- x0 - (f(x0) / dx) # Calculate next value x1
    k[i] <- x1 # Store x1
    # Once the difference between x0 and x1 becomes sufficiently small, output the results.
    if (abs(x1 - x0) < tol) {
      root.approx <- tail(k, n=1)
      res <- list('root_approximation' = root.approx, 'iterations' = k)
      return(res)
    }
    # If Newton-Raphson has not yet reached convergence set x1 as x0 and continue
    x0 <- x1
  }
  print('Too many iterations in method')
}



geracao <- function(n, alpha, beta, p){
  N<-c();G<-c();Gp<-c()
  for (i in 1:n) {
    N[i]<-(rgeom(n = 1, prob = p)+1)
    Gp<-rgamma(n = N[i], shape = alpha, rate = beta)
    G[i]<-sum(Gp)
  }
  NG<-rbind(N,G)
  return(NG)
}



deriv_alpha <- function(alpha){
  n*N_media*log(alpha*(N_media/X_media)) + sum(dat[1,]*log(dat[2,]) - dat[1,]*psigamma(alpha*dat[1,]))
}



estimacao <- function(x){
  n <- dim(x)[2]
  N_media <- mean(x[1,])
  X_media <- mean(x[2,])
  
  assign("n", n, environment(deriv_alpha))
  assign("N_media", N_media, environment(deriv_alpha))
  assign("X_media", X_media, environment(deriv_alpha))
  
  alpha_est <- newton.raphson(
    f = deriv_alpha, a = 0.01, b = 100
  )$root_approximation
  
  # rm("n", envir = environment(deriv_alpha))
  # rm("N_media", envir = environment(deriv_alpha))
  # rm("X_media", envir = environment(deriv_alpha))
  
  beta_est <- alpha_est*(N_media/X_media)
  p_est <- 1/N_media
  est <- c(alpha_est,beta_est,p_est)
  return(est)
}



FVS_Restrita <- function(x){
  
  N <- x[1,]
  X <- x[2,]
  N_media <- mean(N)
  X_media <- mean(X)
  n <- length(N)
  
  beta <- N_media/X_media
  p <- 1/N_media
  
  restrita <- n*N_media*log(beta)-n*N_media+n*log(p)+n*(N_media-1)*log(1-p)+sum(N*log(X)-log(X)-log(gamma(N)))
  result <- list(
    X = X,
    N = N,
    parametros = c(1,beta,p),
    vero_restrita = restrita
  )
  
  return(result)
  
}



FVS_Irrestrita <- function(x){
  
  N <- x[1,]
  X <- x[2,]
  N_media <- mean(N)
  X_media <- mean(X)
  n <- length(N)
  par <- estimacao(dat)
  alpha_est <- par[1]
  beta <- par[2]
  p <- par[3]
  
  irrestrita <- alpha_est*n*N_media*log(beta)-beta*n*X_media+n*log(p)+n*(N_media-1)*log(1-p)+sum(alpha_est*N*log(X)-log(X)-log(gamma(alpha_est*N)))
  
  result <- list(
    X = X
    , N = N
    , parametros = c(alpha_est,beta,p)
    , vero_irrestrita = irrestrita
  )
  
  return(result)
  
}



TRV <- function(x){
  vero_restrita <- FVS_Restrita(x)$vero_restrita
  vero_irrestrita <- FVS_Irrestrita(x)$vero_irrestrita
  TRV_value <- 2*(vero_irrestrita-vero_restrita)
  return(TRV_value)
}



# n = 10 ####

n <- 10
B <- 100
k <- 1000
TRV_value <- trv_boot_mean <- trv_boot <- c()
trv_power <- list()



pb <- progress_bar$new(
  total = k*B#*5
  ,format = "  Processing |:bar| :percent eta: :eta / :current de 1.000.000"
)

for (z in 35:40) {
  for (i in 1:k) {
    set.seed(05062022+i)
    dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    beta_hat <- mean(dat[1,])/mean(dat[2,])
    p_hat <- 1/mean(dat[1,])
    TRV_value[i] <- TRV(dat)
    for (j in 1:B) {
      dat <- geracao(n = n, alpha = 1, beta = beta_hat, p = p_hat)
      trv_boot[j] <- TRV(dat)
      pb$tick()
    }
    trv_boot_mean[i] <- trv_boot %>%
      as.numeric() %>% 
      as_tibble() %>% 
      na.omit() %>% 
      map_dbl(mean)
  }
  trv_power[[z]] <- tibble(
    trv = as.numeric(TRV_value)
    , trv_boot = trv_boot_mean
  ) %>% na.omit() %>% 
    mutate(
      V1 = trv/trv_boot
    ) %>% 
    mutate(
      est = factor(
        ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
        , levels = c('Rejeitou H0', 'Nao Rejeitou H0')
      )
    ) %$% 
    table(est) %>% 
    prop.table()
}



trv_n10 <- trv_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble() %>% 
  bind_cols(parametro = seq(0.1,4, by = 0.1)) %>% 
  dplyr::select(parametro, value)



# n = 20 ####

n <- 20
B <- 100
k <- 1000
TRV_value <- trv_boot_mean <- trv_boot <- c()
trv_power <- list()



pb <- progress_bar$new(
  total = k*B#*5
  ,format = "  Processing |:bar| :percent eta: :eta / :current de 1.000.000"
)

for (z in 35:40) {
  for (i in 1:k) {
    set.seed(05062022+i)
    dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    beta_hat <- mean(dat[1,])/mean(dat[2,])
    p_hat <- 1/mean(dat[1,])
    TRV_value[i] <- TRV(dat)
    for (j in 1:B) {
      dat <- geracao(n = n, alpha = 1, beta = beta_hat, p = p_hat)
      trv_boot[j] <- TRV(dat)
      pb$tick()
    }
    trv_boot_mean[i] <- trv_boot %>%
      as.numeric() %>% 
      as_tibble() %>% 
      na.omit() %>% 
      map_dbl(mean)
  }
  trv_power[[z]] <- tibble(
    trv = as.numeric(TRV_value)
    , trv_boot = trv_boot_mean
  ) %>% na.omit() %>% 
    mutate(
      V1 = trv/trv_boot
    ) %>% 
    mutate(
      est = factor(
        ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
        , levels = c('Rejeitou H0', 'Nao Rejeitou H0')
      )
    ) %$% 
    table(est) %>% 
    prop.table()
}



trv_n20 <- trv_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble() %>% 
  bind_cols(parametro = seq(0.1,4, by = 0.1)) %>% 
  dplyr::select(parametro, value)



# n = 30 ####

n <- 30
B <- 100
k <- 1000
TRV_value <- trv_boot_mean <- trv_boot <- c()
trv_power <- list()



pb <- progress_bar$new(
  total = k*B#*10
  ,format = "  Processing |:bar| :percent eta: :eta / :current de 1.000.000"
)

for (z in 30:40) {
  for (i in 1:k) {
    set.seed(05062022+i)
    dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    beta_hat <- mean(dat[1,])/mean(dat[2,])
    p_hat <- 1/mean(dat[1,])
    TRV_value[i] <- TRV(dat)
    for (j in 1:B) {
      dat <- geracao(n = n, alpha = 1, beta = beta_hat, p = p_hat)
      trv_boot[j] <- TRV(dat)
      pb$tick()
    }
    trv_boot_mean[i] <- trv_boot %>%
      as.numeric() %>% 
      as_tibble() %>% 
      na.omit() %>% 
      map_dbl(mean)
  }
  trv_power[[z]] <- tibble(
    trv = as.numeric(TRV_value)
    , trv_boot = trv_boot_mean
  ) %>% na.omit() %>% 
    mutate(
      V1 = trv/trv_boot
    ) %>% 
    mutate(
      est = factor(
        ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
        , levels = c('Rejeitou H0', 'Nao Rejeitou H0')
      )
    ) %$% 
    table(est) %>% 
    prop.table()
}



trv_n30 <- trv_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble() %>% 
  bind_cols(parametro = seq(0.1,4, by = 0.1)) %>% 
  dplyr::select(parametro, value)



# n = 50 ####

n <- 50
B <- 100
k <- 1000
TRV_value <- trv_boot_mean <- trv_boot <- c()
trv_power <- list()



pb <- progress_bar$new(
  total = k*B*10
  ,format = paste0(
    "  Processing |:bar| :percent eta: :eta / :current de "
    , prettyNum(format(k*B*10, scientific = F), big.mark = '.', decimal.mark = ',')
  )
)



for (z in 30:40) {
  for (i in 1:k) {
    set.seed(05062022+i)
    dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    beta_hat <- mean(dat[1,])/mean(dat[2,])
    p_hat <- 1/mean(dat[1,])
    TRV_value[i] <- TRV(dat)
    for (j in 1:B) {
      dat <- geracao(n = n, alpha = 1, beta = beta_hat, p = p_hat)
      trv_boot[j] <- TRV(dat)
      pb$tick()
    }
    trv_boot_mean[i] <- trv_boot %>%
      as.numeric() %>% 
      as_tibble() %>% 
      na.omit() %>% 
      map_dbl(mean)
  }
  trv_power[[z]] <- tibble(
    trv = as.numeric(TRV_value)
    , trv_boot = trv_boot_mean
  ) %>% na.omit() %>% 
    mutate(
      V1 = trv/trv_boot
    ) %>% 
    mutate(
      est = factor(
        ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
        , levels = c('Rejeitou H0', 'Nao Rejeitou H0')
      )
    ) %$% 
    table(est) %>% 
    prop.table()
}



 trv_n50 <- trv_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble() %>% 
  bind_cols(parametro = seq(0.1,4, by = 0.1)) %>% 
  dplyr::select(parametro, value)


rm(list=ls()[!ls() %in% c('trv_n10','trv_n20','trv_n30','trv_n50')])






bind_cols(trv_n10,trv_n20,trv_n30) %>%
  .[,c(1,2,4,6)] %>% 
  magrittr::set_colnames(c('parametro','n = 10','n = 20','n = 30')) %>% 
  reshape2::melt(id.vars = 'parametro') %>% na.omit() %>% 
  ggplot(aes(x = parametro, y = value, color = variable))+
  geom_line()+geom_point()+
  labs(
    title="Função Poder - Teste TRV"
    , x = "Valor do Parâmetro"
    , y = "Probabilidade de Rejeição"
    , color = 'Tamanho amostral'
  )+
  theme_light()+
  theme(
    legend.position = 'bottom'
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , plot.title = element_text(hjust = 0.5)
  )




