
library(munsell)
library(tidyverse)
library(magrittr)
library(reshape2)
library(progress)

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
  assign("dat", x, environment(deriv_alpha))
  
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



# Teste de Wald

Teste_Wald <- function(x,n){
  wald_par <- estimacao(x)
  alpha <- wald_par[1]
  beta <- wald_par[2]
  p_wald <- wald_par[3]
  res <- NA
  
  err <- 1
  j <- 1
  
  while(err > 1e-16){
    if(j == 1){
      parc <- 5
    } else{
      parc <- res[j-1]
    }
    res[j] <- j^2*psigamma(j*alpha,1)*p_wald*(1-p_wald)^(j-1)
    err <- abs(res[j] - parc)
    j <- j + 1
  }
  
  C <- sum(res)
  # InfoFisher <- matrix(
  #   c(
  #     C, -1/(p_wald*beta), 0
  #     , -1/(p_wald*beta), alpha/(p_wald*beta^2), 0
  #     , 0, 0, 1/(p_wald^2)+1/(p_wald*(1-p_wald))
  #   )
  #   , ncol = 3
  #   ,byrow = T
  # )
  # InfoFisher <- n*InfoFisher
  # diff <- matrix(c(alpha,beta,p_wald) - c(1,N_media/X_media,1/N_media), ncol = 1)
  
  Wald <- as.vector(n*(alpha-1)^2*solve(matrix(C)))
  
  return(Wald)
}



# n = 10 ####

wald_power <- list()
n <- 10
k <- 10000
Wald <- c()
pb <- progress_bar$new(
  total = 50
  ,format = "  Processing |:bar| :percent eta: :eta / :current de :total"
)

for (z in 1:50) {
  for (i in 1:k) {
    wald_dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    Wald[i] <- try(Teste_Wald(wald_dat, n))
  }
  
  wald_power[[z]] <- Wald %>%
    as_tibble() %>%
    mutate(V1 = as.numeric(value)) %>% 
    na.omit() %>% 
    mutate(
      est = factor(
        ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
        , levels = c('Rejeitou H0', 'Nao Rejeitou H0')
      )
    ) %$% 
    table(est) %>% 
    prop.table()
  
  pb$tick()
}



n10 <- wald_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble() %>% 
  bind_cols(parametro = seq(0.1,5, by = 0.1)) %>% 
  dplyr::select(parametro, value)



# n = 20 ####

wald_power <- list()
n <- 20
k <- 10000
Wald <- c()
pb <- progress_bar$new(
  total = 50
  ,format = "  Processing |:bar| :percent eta: :eta / :current de :total"
)

for (z in 1:50) {
  for (i in 1:k) {
    wald_dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    Wald[i] <- try(Teste_Wald(wald_dat, n))
    #pb$tick()
  }
  
  wald_power[[z]] <- Wald %>%
    as_tibble() %>%
    mutate(V1 = as.numeric(value)) %>% 
    na.omit() %>% 
    mutate(
      est = factor(
        ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
        , levels = c('Rejeitou H0', 'Nao Rejeitou H0')
      )
    ) %$% 
    table(est) %>% 
    prop.table()
  
  pb$tick()
}



n20 <- wald_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble()



# n = 30 ####

wald_power <- list()
n <- 30
k <- 10000
Wald <- c()
pb <- progress_bar$new(
  total = 50
  ,format = "  Processing |:bar| :percent eta: :eta / :current de :total"
)

for (z in 1:50) {
  for (i in 1:k) {
    wald_dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    Wald[i] <- try(Teste_Wald(wald_dat, n))
    #pb$tick()
  }
  
  wald_power[[z]] <- Wald %>%
    as_tibble() %>%
    mutate(V1 = as.numeric(value)) %>% 
    na.omit() %>% 
    mutate(
      est = factor(
        ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
        , levels = c('Rejeitou H0', 'Nao Rejeitou H0')
      )
    ) %$% 
    table(est) %>% 
    prop.table()
  
  pb$tick()
}



n30 <- wald_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble()



# n = 50 ####

wald_power <- list()
n <- 50
k <- 10000
Wald <- c()
pb <- progress_bar$new(
  total = 50
  ,format = "  Processing |:bar| :percent eta: :eta / :current de :total"
)

for (z in 1:50) {
  for (i in 1:k) {
    wald_dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    Wald[i] <- try(Teste_Wald(wald_dat, n))
    #pb$tick()
  }
  
  wald_power[[z]] <- Wald %>%
    as_tibble() %>%
    mutate(V1 = as.numeric(value)) %>% 
    na.omit() %>% 
    mutate(
      est = factor(
        ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
        , levels = c('Rejeitou H0', 'Nao Rejeitou H0')
      )
    ) %$% 
    table(est) %>% 
    prop.table()
  
  pb$tick()
}



n50 <- wald_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble()



wald_power_sb <- bind_cols(n10,n20,n30,n50) %>% 
  magrittr::set_colnames(c('parametro','n = 10','n = 20','n = 30','n = 50')) %>% 
  reshape2::melt(id.vars = 'parametro') %>% 
  ggplot(aes(x = parametro, y = value, color = variable))+
  geom_line()+geom_point()+
  labs(
    title="Fun????o Poder - Teste Wald"
    , x = "Valor do Par??metro"
    , y = "Probabilidade de Rejei????o"
    , color = 'Tamanho amostral'
  )+
  theme_light()+
  theme(
    legend.position = 'bottom'
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , plot.title = element_text(hjust = 0.5)
  )
wald_power_sb

rm(list = ls()[ls() != 'wald_power_sb'])


