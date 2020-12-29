# Handmade EM4MoG ---------------------------------------------------------

vnorm <- function(y, p, mu, sigma) p * dnorm(y, mu, sigma)
vect_wnorm <- Vectorize(vnorm, c("p", "mu", "sigma"))
vect_norm = Vectorize(dnorm, c("mu", "sigma"))

handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = T)
{
  # vectorized likelihood
  like     <- apply(vect_wnorm(y, p, mu, sigma), 1, sum)
  deviance <- -2 * sum(log(like))
  
  # sampling plot colors
  color_list <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  n    <- length(p)
  cols <- sample(color_list, n)

  # creating matrix of iteration values
  n_col    <- 2 + 3 * n
  res      <- matrix(NA, n_iter + 1, n_col)
  res[1,]  <- c(0, p, mu, sigma, deviance)
  
  for (iter in 1:n_iter) {
    # E step
    d   <- vect_wnorm(y, p, mu, sigma) # matrix of responsabilities
    den <- apply(d, 1, sum) # normalizing factor
    r   <- apply(d, 2, function(x) x / den) # normalizing each column
    
    # M step
    p <- apply(r, 2, mean)
    mu <- apply(r, 2, function(x) sum(x * y) / sum(x))
    sigma <- apply(r, 2, function(x) sum(x * (y^2)) / sum(x))
    sigma <- sqrt(sigma - mu^2)
    
    # -2 x log-likelihood (a.k.a. deviance)
    like     <- apply(vect_wnorm(y, p, mu, sigma), 1, sum)
    deviance <- -2 * sum(log(like))
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
    # Plot
    if (plot_flag){
      hist(y, prob = T, breaks = 30, col = gray(.8), border = NA, 
           main = "", xlab = paste("EM Iteration: ", iter, "/", n_iter, sep = ""))
      set.seed(123)
      # compute max component for each obs
      norms <- vect_norm(y, mu, sigma)
      cols_idxs <- apply(norms, 1, which.max)
      curve_fun <- function(x) sum(vect_wnorm(x, p, mu, sigma))
      
      points(jitter(y), rep(0,length(y)), 
             pch = 19, cex = .6, 
             col = cols[cols_idxs])
      curve(sapply(x, curve_fun),
            lwd = 4, col = rgb(0,0,0,.5), add = TRUE)
      Sys.sleep(1.5)
    }
  }
  res <- data.frame(res)
  # function to generate strings
  append_str <- function(str) sapply(1:n, function(x) paste(str, x, sep=''))
  p_names <- append_str("p")
  mu_names <- append_str("mu")
  sigma_names <- append_str("sigma")
  names(res) <- c("iteration", p_names, mu_names, sigma_names, "deviance")
  out <- list(p = p, mu = mu, sigma = sigma, deviance = deviance, res = res)
  return(out)
}

# Bart Simpson density
bart_simpson <- function(x) 0.5*dnorm(x, 0, 1) +
                            0.1*dnorm(x, -1.0, 0.1) + 0.1*dnorm(x, -0.5, 0.1) +
                            0.1*dnorm(x, 0.0, 0.1) + 0.1*dnorm(x, 0.5, 0.1) +
                            0.1*dnorm(x, 1.0, .1)

# create a random sample
sample_sim <- function(n) rnormmix(n,
                                   lambda = c(0.5, rep(0.1,5)),
                                   mu     = c(0, ((0:4)/2)-1),
                                   sigma  = c(1, rep(0.1,5)))

n_1 = 500 # non-asymptotic
n_2 = 100000 # asymptotic
n_breaks = 100

# hist(sample_sim(n_1), prob=T, breaks=n_breaks) # uncomment to plot
hist(sample_sim(n_2), prob=T, breaks=n_breaks)

XX <- sample_sim(n_2)
hm.fit <- handmade.em(XX, p=rep(1, 5)/5,
                      mu=runif(5, 1, 4),
                      sigma=runif(5, 0, 2),
                      n_iter=10)

log_likelihood <- function(x, hm.fit) {
  like <- apply(vect_wnorm(x, hm.fit$p, hm.fit$mu, hm.fit$sigma), 1, sum)
  log_like <- sum(log(like)) # log likelihood
  return(log_like)
}


init_params <- function(k) {
  # takes the data and the number of components and returns random
  # initial parameters, data are used to compute min/max to
  # generate mu and sigma
  p <- runif(k, 0, 1)
  p <- p / sum(p) # normalizing
  mu <- runif(k, -1, 1)
  sigma <- runif(k, 0, 1)
  return(list(p = p, mu = mu, sigma = sigma))
}

sample_split <- function(x, ratio) {
  # x : data
  # ratio : float (represents length(train)/)
  # returns a random split of x in train and test set based on ratio
  n <- length(x)
  n_train <- as.integer(n * ratio)
  train_idx <- sample(1:n, size = n_train) # sample indexes for trainset
  test_idx <- setdiff(1:n, train_idx)
  return(list(train = x[train_idx], test = x[test_idx]))
}

kfold_to_split <- function(kfold, k) {
  # kfold : is the result of createFolds() with names Fold1, Fold2, ...
  # returns a list with names train and test, setting Fold_k as test set
  # and the rest as training set
  idxs <- setdiff(1:length(kfold), k)
  train_set <- unlist(kfold[idxs], use.names = F)
  return(list(train = train_set, test = kfold[[k]]))
}

kfold_cv <- function(x, k, em_run) {
  # x : data
  # k : number of folds
  # returns the average risk running a kfold cross validation on x
  kfold <- createFolds(x, k) # named list of Fold1, Fold2, ...
  outs <- rep(NA, k)
  
  for (i in 1:k) {
    x_split <- kfold_to_split(kfold, i) # named list of train/test
    hm.fit <- em_run(x_split$train)
    log_like <- 1/n * log_likelihood(x_split$test, hm.fit)
    outs[i] <- log_like
  }
  
  return(mean(outs))
}


EM_simulation <- function(n, M, k_max, n_iter = 20) {
  # n : sample size
  # M : number of simulations
  # k_max : maximum number of components
  # n_iter : number of iterations of EM
  
  results <- c()
  
  for (n_sim in 1:M) {
    X <- sample_sim(n) # simulating the Bart

    # model selection evaluations
    model_eval <- matrix(nrow = 8, ncol = k_max)
    
    for (k in 1:k_max) {
      params <- init_params(k) # initial random parameters
      
      # hem.fit instantiated with fixed parameters to make multiple calls within the loop
      em_run <- function(x) handmade.em(X, params$p, params$mu, params$sigma, n_iter, plot = F)
      
      # first run EM algorithm on the whole dataset to compute AIC and BIC
      hm.fit <- em_run(X)
      
      # AIC and BIC evaluation
      log_like <- log_likelihood(X, hm.fit)
      model_eval[1,k] <- 2 * (k - log_like)
      model_eval[2,k] <- log(n) * k - 2 * log_like
      
      # cross validation on different split ratios
      ratios = c(.5, .7, .3)
      
      for (i in 1:length(ratios)) {
        x_split <- sample_split(x, ratios[i])
        hm.fit <- em_run(x_split$train)
        log_like <- 1/n * log_likelihood(x_split$test, hm.fit)
        model_eval[i+2,k] <- log_like
      }
      
      # 5-fold cross validation
      model_eval[6,k] <- kfold_cv(x, 5, em_run)
      
      # 10-fold cross validation
      model_eval[7,k] <- kfold_cv(x, 5, em_run)
    }
  
    cat("sim =", n_sim, "\n")
    print(model_eval)
  }
}
