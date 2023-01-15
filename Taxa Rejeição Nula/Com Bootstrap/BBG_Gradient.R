
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

Teste_Gradient <- function(x,n){
  
  N <- x[1,]
  X <- x[2,]
  
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



# Alpha = 0.9 ####



# n = 10

n <- 10
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.9_10 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.9,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 20

n <- 20
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.9_20 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.9,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 30

n <- 30
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.9_30 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.9,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 40

n <- 40
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.9_40 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.9,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 50

n <- 50
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.9_50 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.9,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()

score_0.9_10
score_0.9_20
score_0.9_30
score_0.9_50


tab_result_0.9 <- rbind(
  score_0.9_10 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.9_20 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.9_30 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.9_40 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.9_50 %>% as.matrix.noquote() %>% unname() %>% t()
)
row.names(tab_result_0.9) <- c('n = 10', 'n = 20', 'n = 30', 'n = 40', 'n = 50')
colnames(tab_result_0.9) <- c('Rejeitou H0', 'Nao Rejeitou H0')
tab_result_0.9



# Alpha = 0.95 ####



# n = 10

n <- 10
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.95_10 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 20

n <- 20
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.95_20 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 30

n <- 30
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.95_30 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 40

n <- 40
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.95_40 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 50

n <- 50
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.95_50 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.95,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()

score_0.95_10
score_0.95_20
score_0.95_30
score_0.95_50


tab_result_0.95 <- rbind(
  score_0.95_10 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.95_20 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.95_30 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.95_40 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.95_50 %>% as.matrix.noquote() %>% unname() %>% t()
)
row.names(tab_result_0.95) <- c('n = 10', 'n = 20', 'n = 30', 'n = 40', 'n = 50')
colnames(tab_result_0.95) <- c('Rejeitou H0', 'Nao Rejeitou H0')
tab_result_0.95




# Alpha = 0.99 ####



# n = 10

n <- 10
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.99_10 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.99,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 20

n <- 20
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.99_20 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.99,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 30

n <- 30
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.99_30 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.99,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 40

n <- 40
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.99_40 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.99,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()



# n = 50

n <- 50
k <- 10000
grad <- c()
pb <- progress_bar$new(total = k)

for (i in 1:k) {
  set.seed(28092022+i)
  grad_dat <- geracao(n = n, alpha = 1, beta = 2, p = 0.5)
  grad[i] <- try(Teste_Gradient(grad_dat, n))
  pb$tick()
}

score_0.99_50 <- grad %>%
  as_tibble() %>%
  mutate(V1 = as.numeric(value)) %>% 
  na.omit() %>% 
  mutate(
    est = ifelse(V1 > qchisq(0.99,1), 'Rejeitou H0', 'Nao Rejeitou H0')
  ) %$% 
  table(est) %>% 
  prop.table()

score_0.99_10
score_0.99_20
score_0.99_30
score_0.99_50


tab_result_0.99 <- rbind(
  score_0.99_10 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.99_20 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.99_30 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.99_40 %>% as.matrix.noquote() %>% unname() %>% t()
  , score_0.99_50 %>% as.matrix.noquote() %>% unname() %>% t()
)
row.names(tab_result_0.99) <- c('n = 10', 'n = 20', 'n = 30', 'n = 40', 'n = 50')
colnames(tab_result_0.99) <- c('Rejeitou H0', 'Nao Rejeitou H0')
tab_result_0.99

list(
  Alpha_0.9 = tab_result_0.9
  , Alpha_0.95 = tab_result_0.95
  , Alpha_0.99 = tab_result_0.99
)


