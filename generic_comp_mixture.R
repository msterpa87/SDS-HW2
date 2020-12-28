M   <- 100
# M   <- 1000
out <- rep(NA, M)
for (m in 1:M){
  lab = sample(c("Female", "Male"), size = 1, prob = c(.4,.6))
  if (lab == "Female") 
    out[m] <- rnorm(n = 1, mean = 162, sd = sqrt(5))
  else 
    out[m] <- rnorm(n = 1, mean = 176, sd = sqrt(8))
}
hist(out, main = "", xlab = "")

# Handmade EM4MoG ---------------------------------------------------------

vnorm <- function(p, mu, sigma, y) p * dnorm(y, mu, sigma)
vect_wnorm <- Vectorize(vnorm, c("p", "mu", "sigma"))

handmade.em <- function(y, p, mu, sigma, n_iter, plot_flag = T)
{
  # vectorized likelihood
  like     <- vect_wnorm(y, p, mu, sigma)
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
    p[1]     <- mean(r1)
    mu[1]    <- sum(r1*y)/sum(r1)
    sigma[1] <-sqrt( sum(r1*(y^2))/sum(r1) - (mu[1])^2 )
    p[2]     <- 1 - p[1]
    mu[2]    <- sum((r2)*y)/sum((r2))
    sigma[2] <- sqrt(sum(r2*(y^2))/sum(r2) - (mu[2])^2)
    
    # -2 x log-likelihood (a.k.a. deviance)
    like     <- p[1]*dnorm(y, mu[1], sigma[1]) + p[2]*dnorm(y, mu[2], sigma[2])
    deviance <- -2*sum( log(like) )
    
    # Save
    res[iter+1,] <- c(iter, p, mu, sigma, deviance)
    
    # Plot
    if (plot_flag){
      hist(y, prob = T, breaks = 30, col = gray(.8), border = NA, 
           main = "", xlab = paste("EM Iteration: ", iter, "/", n_iter, sep = ""))
      set.seed(123)
      points(jitter(y), rep(0,length(y)), 
             pch = 19, cex = .6, 
             col = cols[ (dnorm(y,mu[1],sigma[1]) > dnorm(y,mu[2],sigma[2])) + 1])
      curve(p[1]*dnorm(x,mu[1],sigma[1]) + p[2]*dnorm(x,mu[2],sigma[2]),
            lwd = 4, col = rgb(0,0,0,.5), add = TRUE)
      Sys.sleep(1.5)
    }
  }
  res <- data.frame(res)
  names(res) <- c("iteration","p1","p2","mu1","mu2","sigma1","sigma2","deviance")
  out <- list(parameters = c(p = p, mu = mu, sigma = sigma), deviance = deviance, res = res)
  return(out)
}