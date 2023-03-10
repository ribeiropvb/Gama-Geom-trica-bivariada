

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


data <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)



Teste_Gradient <- function(x){
  
  N <- x[1,]
  X <- x[2,]
  n <- length(X)
  
  pars <- estimacao(x)
  
  alpha <- pars[1]
  beta <- pars[2]
  p_grad <- pars[3]
  
  err <- 1
  j <- 1
  res <- NA
  
  while(err > 1e-16){
    if(j == 1){
      parc <- 5
    } else{
      parc <- res[j-1]
    }
    res[j] <- j^2*psigamma(j*alpha,1)*p_grad*(1-p_grad)^(j-1)
    err <- abs(res[j] - parc)
    j <- j + 1
  }
  
  C <- sum(res)
  InfoFisher <- matrix(
    c(
      C, -1/(p_grad*beta), 0
      , -1/(p_grad*beta), alpha/(p_grad*beta^2), 0
      , 0, 0, 1/(p_grad^2)+1/(p_grad*(1-p_grad))
    )
    , ncol = 3
    ,byrow = T
  )
  InfoFisher <- n*InfoFisher
  
  score <- n*mean(N)*log(mean(N)/mean(X))+sum(N*log(X)-N*psigamma(N,0))
  
  res <- score*(alpha-1)
  
  return(res)
}



# n = 10 ####

n <- 10
B <- 100
k <- 1000
grad_value <- grad_boot_mean <- grad_boot <- c()
grad_power <- list()



pb <- progress_bar$new(
  total = k*B*20
  ,format = paste0(
    "  Processing |:bar| :percent eta: :eta / :current de "
    , prettyNum(format(k*B*10, scientific = F), big.mark = '.', decimal.mark = ',')
  )
)

for (z in 20:40) {
  for (i in 1:k) {
    set.seed(05062022+i)
    dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    beta_hat <- mean(dat[1,])/mean(dat[2,])
    p_hat <- 1/mean(dat[1,])
    grad_value[i] <- try(Teste_Gradient(dat), T)
    for (j in 1:B) {
      dat <- geracao(n = n, alpha = 1, beta = beta_hat, p = p_hat)
      grad_boot[j] <- try(Teste_Gradient(dat), T)
      pb$tick()
    }
    grad_boot_mean[i] <- grad_boot %>%
      as.numeric() %>% 
      as_tibble() %>% 
      na.omit() %>% 
      map_dbl(mean)
  }
  grad_power[[z]] <- tibble(
    trv = as.numeric(grad_value)
    , grad_boot = grad_boot_mean
  ) %>% na.omit() %>% 
    mutate(
      V1 = trv/grad_boot
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



grad_n10 <- grad_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble() %>% 
  bind_cols(parametro = seq(0.1,4, by = 0.1)) %>% 
  dplyr::select(parametro, value)



# n = 20 ####

n <- 20
B <- 100
k <- 1000
grad_value <- grad_boot_mean <- grad_boot <- c()
grad_power <- list()



pb <- progress_bar$new(
  total = k*B*20
  ,format = paste0(
    "  Processing |:bar| :percent eta: :eta / :current de "
    , prettyNum(format(k*B*20, scientific = F), big.mark = '.', decimal.mark = ',')
  )
)

for (z in 20:40) {
  for (i in 1:k) {
    set.seed(05062022+i)
    dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    beta_hat <- mean(dat[1,])/mean(dat[2,])
    p_hat <- 1/mean(dat[1,])
    grad_value[i] <- try(Teste_Gradient(dat), T)
    for (j in 1:B) {
      dat <- geracao(n = n, alpha = 1, beta = beta_hat, p = p_hat)
      grad_boot[j] <- try(Teste_Gradient(dat), T)
      pb$tick()
    }
    grad_boot_mean[i] <- grad_boot %>%
      as.numeric() %>% 
      as_tibble() %>% 
      na.omit() %>% 
      map_dbl(mean)
  }
  grad_power[[z]] <- tibble(
    trv = as.numeric(grad_value)
    , grad_boot = grad_boot_mean
  ) %>% na.omit() %>% 
    mutate(
      V1 = trv/grad_boot
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



grad_n20 <- grad_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble() %>% 
  bind_cols(parametro = seq(0.1,4, by = 0.1)) %>% 
  dplyr::select(parametro, value)



# n = 30 ####

n <- 30
B <- 100
k <- 1000
grad_value <- grad_boot_mean <- grad_boot <- c()
grad_power <- list()



pb <- progress_bar$new(
  total = k*B*20
  ,format = paste0(
    "  Processing |:bar| :percent eta: :eta / :current de "
    , prettyNum(format(k*B*20, scientific = F), big.mark = '.', decimal.mark = ',')
  )
)

for (z in 20:40) {
  for (i in 1:k) {
    set.seed(05062022+i)
    dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    beta_hat <- mean(dat[1,])/mean(dat[2,])
    p_hat <- 1/mean(dat[1,])
    grad_value[i] <- try(Teste_Gradient(dat), T)
    for (j in 1:B) {
      dat <- geracao(n = n, alpha = 1, beta = beta_hat, p = p_hat)
      grad_boot[j] <- try(Teste_Gradient(dat), T)
      pb$tick()
    }
    grad_boot_mean[i] <- grad_boot %>%
      as.numeric() %>% 
      as_tibble() %>% 
      na.omit() %>% 
      map_dbl(mean)
  }
  grad_power[[z]] <- tibble(
    trv = as.numeric(grad_value)
    , grad_boot = grad_boot_mean
  ) %>% na.omit() %>% 
    mutate(
      V1 = trv/grad_boot
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



grad_n30 <- grad_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble() %>% 
  bind_cols(parametro = seq(0.1,4, by = 0.1)) %>% 
  dplyr::select(parametro, value)



# n = 50 ####

n <- 50
B <- 100
k <- 1000
grad_value <- grad_boot_mean <- grad_boot <- c()
grad_power <- list()



pb <- progress_bar$new(
  total = k*B*20
  ,format = paste0(
    "  Processing |:bar| :percent eta: :eta / :current de "
    , prettyNum(format(k*B*20, scientific = F), big.mark = '.', decimal.mark = ',')
  )
)

for (z in 20:40) {
  for (i in 1:k) {
    set.seed(05062022+i)
    dat <- geracao(n = n, alpha = z/10, beta = 2, p = 0.5)
    beta_hat <- mean(dat[1,])/mean(dat[2,])
    p_hat <- 1/mean(dat[1,])
    grad_value[i] <- try(Teste_Gradient(dat), T)
    for (j in 1:B) {
      dat <- geracao(n = n, alpha = 1, beta = beta_hat, p = p_hat)
      grad_boot[j] <- try(Teste_Gradient(dat), T)
      pb$tick()
    }
    grad_boot_mean[i] <- grad_boot %>%
      as.numeric() %>% 
      as_tibble() %>% 
      na.omit() %>% 
      map_dbl(mean)
  }
  grad_power[[z]] <- tibble(
    trv = as.numeric(grad_value)
    , grad_boot = grad_boot_mean
  ) %>% na.omit() %>% 
    mutate(
      V1 = trv/grad_boot
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



grad_n50 <- grad_power %>% 
  map(function(x)  x %>% as.matrix.noquote() %>% unname() %>% t() %>% round(4) %>% .[,1]) %>% 
  unlist() %>% as_tibble() %>% 
  bind_cols(parametro = seq(0.1,4, by = 0.1)) %>% 
  dplyr::select(parametro, value)



bind_cols(
  grad_n10
  , grad_n20 %>% .$value
  , grad_n30 %>% .$value
  , grad_n50 %>% .$value
) %>% 
  magrittr::set_colnames(
    c(
      'parametro'
      ,'n = 10'
      ,'n = 20'
      ,'n = 30'
      ,'n = 50'
    )
  ) %>% 
  reshape2::melt(id.vars = 'parametro') %>% 
  ggplot(aes(x = parametro, y = value, color = variable))+
  geom_line()+geom_point()+
  labs(
    title="Fun????o Poder - Teste Score"
    , x = "Probabilidade de Rejei????o"
    , y = "Valor do Par??metro"
    , color = 'Tamanho amostral'
  )+
  theme_light()+
  theme(
    legend.position = 'bottom'
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , plot.title = element_text(hjust = 0.5)
  )


rm(list=ls()[!ls() %in% c('grad_n10','grad_n20','grad_n30','grad_n50')])


