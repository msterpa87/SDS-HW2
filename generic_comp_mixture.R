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
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma), deviance = deviance, res = res)
  return(out)
}

# Bart Simpson density
bart_simpson <- function(x) 0.5*dnorm(x, 0, 1) +
                            0.1*dnorm(x, -1.0, 0.1) + 0.1*dnorm(x, -0.5, 0.1) +
                            0.1*dnorm(x, 0.0, 0.1) + 0.1*dnorm(x, 0.5, 0.1) +
                            0.1*dnorm(x, 1.0, .1)

sample_sim <- function(n) rnormmix(n,
                                   lambda = c(0.5, rep(0.1,5)),
                                   mu     = c(0, ((0:4)/2)-1),
                                   sigma  = c(1, rep(0.1,5)))

#curve(bart_simpson, col=rgb(1,0,0,0,4), lwd=3, n=500)

# n_1: hist(sample_sim(500), prob=T, breaks=100)
# n_2: hist(sample_sim(100000), prob=T, breaks=100)
