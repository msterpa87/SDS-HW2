# Handmade EM4MoG ---------------------------------------------------------

library(mixtools)
library(caret)
library(KScorrect)

vnorm <- function(y, p, mu, sigma) p * dnorm(y, mu, sigma)
vect_wnorm <- Vectorize(vnorm, c("p", "mu", "sigma"))
log_vnorm <- function(y, p, mu, sigma) p * dnorm(y, mu, sigma, log = TRUE)
vect_log_wnorm <- Vectorize(vnorm, c("p", "mu", "sigma"))
vect_norm = Vectorize(dnorm, c("mean", "sd"))

plot_fit <- function(y, p, mu, sigma, xlab = "", breaks = 50) {
  # sampling plot colors
  color_list <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  k    <- length(p)
  cols <- sample(color_list, k)
  
  hist(y, prob = T, breaks = breaks, col = gray(.8), border = NA, 
       main = "", xlab = xlab)
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
}

handmade.em <- function(y, p, mu, sigma, breaks = 50, n_iter = 1000,
                        threshold = 1e-3, plot = F, plot_freq = 25)
{
  # vectorized likelihood
  like     <- apply(vect_log_wnorm(y, p, mu, sigma), 1, sum)
  deviance <- -2 * sum(like)
  k <- length(p)
  
  # creating matrix of iteration values
  n_col    <- 2 + 3 * k
  res      <- matrix(NA, n_iter + 1, n_col)
  res[1,]  <- c(0, p, mu, sigma, deviance)
  
  for (iter in 1:n_iter) {
    # E step
    d   <- vect_wnorm(y, p, mu, sigma) # matrix of responsabilities
    den <- apply(d, 1, sum) # normalizing factor
    r <- d / den # normalizing each column
    r[is.nan(r)] = 0

    # M step
    p <- apply(r, 2, mean)
    mu <- apply(r, 2, function(x) sum(x * y) / sum(x))
    #for (i in 1:k) {
    #  mu[i] <- sum(r[,i] * y) / sum(r[,i])
    #}

    for (i in 1:k) sigma[i] <- sqrt(sum(r[,i] * (y - mu[i])^2) / sum(r[,i]))
    #sigma <- apply(1:k, 2, function(i) sqrt(sum(r[,i] * (y - mu[i])^2) / sum(r[,i]) ) )
    #sigma <- sqrt(sigma - mu^2)
    
    # -2 x log-likelihood (a.k.a. deviance)
    new_like <- apply(vect_log_wnorm(y, p, mu, sigma), 1, sum)
    
    # stopping criteria
    diff <- abs(sum(like) - sum(new_like))
    if (!is.na(diff) && diff < threshold) break
    like     <- new_like
    deviance <- -2 * sum(like)
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
    # Plot
    if (plot && iter %% plot_freq == 0) {
      xlab <- paste("EM Iteration: ", iter, "/", n_iter, sep = "")
      plot_fit(y, p, mu, sigma, xlab)
      Sys.sleep(1)
    }
  }
  res <- data.frame(res)
  
  # function to generate strings
  append_str <- function(str) sapply(1:k, function(x) paste(str, x, sep=''))
  p_names <- append_str("p")
  mu_names <- append_str("mu")
  sigma_names <- append_str("sigma")
  names(res) <- c("iteration", p_names, mu_names, sigma_names, "deviance")
  out <- list(p = p, mu = mu, sigma = sigma, deviance = deviance, res = res, iter = iter)
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

# hist(sample_sim(n_1), prob=T, breaks=n_breaks) # uncomment to plot
#hist(sample_sim(n_2), prob=T, breaks=n_breaks)

log_likelihood <- function(x, hm.fit) {
  # x : vector
  # hm.fit : output of handmade.em function
  like <- apply(vect_log_wnorm(x, hm.fit$p, hm.fit$mu, hm.fit$sigma), 1, sum)
  log_like <- sum(like)
  return(log_like)
}

init_kmeans <- function(y, k) {
  # takes the data and the number of components and returns
  # initial parameters (proportions, mean, std) running the
  # kmeans algotithm
  labels <- kmeans(y, k)$cluster
  mu <- unlist(lapply(1:k, function(i) mean(y[labels == i])))
  sigma <- unlist(lapply(1:k, function(i) sd(y[labels == i])))
  sigma[is.na(sigma)] = 0
  p <- unlist(lapply(1:k, function(i) sum(labels == i) / length(labels)))
  return(list(p = p, mu = mu, sigma = sigma))
}

init_params <- function(k) {
  # takes the data and the number of components and returns random
  # initial parameters, data are used to compute min/max to
  # generate mu and sigma
  p <- runif(k, 0, 1)
  p <- p / sum(p)
  mu <- runif(k, 0, 1)
  sigma <- runif(k, 0, 1)
  return(list(p = p, mu = mu, sigma = sigma))
}

init_params_groups <- function(y, k) {
  # takes the data and the number of components and returns random
  # initial parameters, data are used to compute min/max to
  # generate mu and sigma
  folds <- createFolds(y, k) # group of indexes
  folds <- rapply(folds, function(x) y[x], how = "replace") # groups of values
  p <- rep(1 / k, k) # same proportion for each component
  p <- p / sum(p) # normalizing
  mu <- as.vector(rapply(folds, mean)) # mean of each random group
  sigma <- as.vector(rapply(folds, sd))
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

kfold_to_split <- function(XX, kfold, k) {
  # kfold : is the result of createFolds() with names Fold1, Fold2, ...
  # returns a list with names train and test, setting Fold_k as test set
  # and the rest as training set
  train_fold_idxs <- setdiff(1:length(kfold), k)
  train_idxs <- unlist(kfold[train_fold_idxs], use.names = F)
  test_idxs <- kfold[[k]]
  return(list(train = XX[train_idxs], test = XX[test_idxs]))
}

kfold_cv <- function(XX, n_fold, k) {
  # x : data
  # k : number of folds
  # returns the average risk running a kfold cross validation on x
  kfold <- createFolds(XX, n_fold) # named list of Fold1, Fold2, ...
  outs <- rep(NA, n_fold)
  
  for (i in 1:n_fold) {
    x_split <- kfold_to_split(XX, kfold, i) # named list of train/test
    params <- init_kmeans(x_split$train, k)
    hm.fit <- handmade.em(x_split$train, params$p, params$mu, params$sigma)
    outs[i] <- log_likelihood(x_split$test, hm.fit) / length(x_split$test)
  }
  return(mean(outs[!is.nan(outs)]))
}

print_params <- function(p, mu, sigma) cat("\np =", p, "\nmu =", mu, "\nsigma =", sigma, "\n")

EM_simulation <- function(n, M, k_max, subdivisions = 200, verbose = TRUE) {
  # n : sample size
  # M : number of simulations
  # k_max : maximum number of components
  # n_iter : number of iterations of EM
  
  result <- matrix(nrow = 8, ncol = M)
  model_performance <- array(dim = c(8, k_max, M))

  for (n_sim in 1:M) {
    
    # sample from Bart
    XX <- sample_sim(n)
    
    # model evaluation matrix
    model_eval <- matrix(nrow = 8, ncol = k_max)

    for (k in 1:k_max) {
      if(verbose) cat("[sim=", n_sim, ", k=", k, "]\n", sep='')
      
      # first run EM algorithm on the whole dataset to compute AIC and BIC
      params <- init_kmeans(XX, k) # initial random parameters
      hm.fit <- handmade.em(XX, params$p, params$mu, params$sigma)
      
      # AIC and BIC evaluation
      log_like <- log_likelihood(XX, hm.fit)
      model_eval[1,k] <- 2 * log_like - 2 * k * 3
      model_eval[2,k] <- log_like - log(n)/n * k * 3
      
      # cross validation on different split ratios
      ratios = c(.5, .7, .3)
      
      for (i in 1:3) {
        x_split <- sample_split(XX, ratios[i])
        params <- init_kmeans(x_split$train, k) # initial parameters
        hm.fit <- handmade.em(x_split$train, params$p, params$mu, params$sigma)
        log_like <- log_likelihood(x_split$test, hm.fit) / length(x_split$test)
        model_eval[i+2,k] <- log_like
      }
      
      # 5-fold cross validation
      model_eval[6,k] <- kfold_cv(XX, 5, k)
      
      # 10-fold cross validation
      model_eval[7,k] <- kfold_cv(XX, 10, k)
      
      # Wasserstein score
      x_split <- sample_split(XX, .5) # split data in half
      params <- init_kmeans(x_split$train, k)
      hm.fit <- handmade.em(x_split$train, params$p, params$mu, params$sigma) # get MLE from trainset
      f_k <- function(p) qmixnorm(p, mean = hm.fit$mu,
                                  sd = hm.fit$sigma,
                                  pro = hm.fit$p,
                                  expand = .2) # quantile function on the MLE
      test_ecdf <- ecdf(x_split$test) # ecdf of the test set
      f_test <- function(p) as.numeric(quantile(test_ecdf, probs = p)) # quantile function of testset ecdf
      f_abs <- function(p) abs(f_k(p) - f_test(p)) # integrand
      model_eval[8,k] <- integrate(f_abs, lower = 0, upper = 1,
                                   subdivisions = subdivisions,
                                   stop.on.error = FALSE)$value
    }

    # take best k for each mode evaluation and save it in the result matrix
    model_eval[is.nan(model_eval)] = -Inf # eliminating possible NaN values
    model_performance[,,n_sim] = model_eval
    result[,n_sim] <- tryCatch(unlist(apply(model_eval, 1, which.max)))  # select the best k for each simulation
  }
  
  return(list("best_k" = result, "performance" = model_performance))
}

accuracy <- function(k_matrix, k) {
  result <- unlist(apply(k_matrix, 1, function(x) sum(x == k) / length(x)))
  #names(result) <- c("AIC", "BIC", "50/50-SS", "70/30-SS", "30/70-SS", "5-Fold-CV", "10-Fold-CV", "Wasserstein")
  return(result)
}

make.accuracy.table <- function(k_matrix, k_max) {
  model.names <- c("AIC", "BIC", "50/50-SS", "70/30-SS", "30/70-SS", "5-Fold-CV", "10-Fold-CV", "Wasserstein")
  accuracy.table <- ldply(1:k_max, function(i) accuracy(k_matrix, i))
  colnames(accuracy.table) <- model.names
  k <- 1:8
  return(cbind(k, accuracy.table))
}

n_1 <- 100
n_2 <- 1000
M <- 5
k_max <- 8
large.result <- EM_simulation(n_1, M, k_max)
small.result <- EM_simulation(n_2, M, k_max)

make.accuracy.table(large.result$best_k, k_max) %>% gt() %>%
  tab_header(
    title = "Model Accuracy (n = 100)"
  ) %>%
  fmt_percent(
    columns = 2:9,
    decimals = 0
  )