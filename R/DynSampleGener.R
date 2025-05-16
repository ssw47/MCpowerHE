# Load necessary library
library(copula)

Generating_Sample <- function(
    endpoints,
    copula_type = "Clayton",    # Options: "Clayton", "Frank", "Gumbel", "Gaussian"
    copula_param = NULL,        # For Archimedes ("Clayton", "Frank", "Gumbel"): Input Kendall's tau; For Gaussian: correlation matrix
    Fellow_up.Time = 180,
    N.Super = 20000
){
  N.endpoints <- length(endpoints)
  cop_dim <- N.endpoints

  # Create the copula
  if (copula_type %in% c("Clayton", "Frank", "Gumbel")){
    # For Archimedes copulas, expect a single parameter (Kendall's tau)
    if (is.null(copula_param)){
      stop("Please provide a copula_param (Kendall's tau) for the Archimedean copula.")
    }
    # Convert Kendall's tau to copula parameter theta
    tau <- copula_param
    if (copula_type == "Clayton"){
      theta <- 2 * tau / (1 - tau)
      cop <- claytonCopula(param = theta, dim = cop_dim)
    } else if (copula_type == "Frank"){
      # Use numerical inversion for theta
      cop_temp <- frankCopula(dim = cop_dim)
      theta <- iTau(cop_temp, tau)
      cop <- frankCopula(param = theta, dim = cop_dim)
    } else if (copula_type == "Gumbel"){
      theta <- 1 / (1 - tau)
      cop <- gumbelCopula(param = theta, dim = cop_dim)
    }
  } else if (copula_type == "Gaussian"){
    # For Gaussian copula, expect a correlation matrix
    if (is.null(copula_param)){
      stop("Please provide a copula_param (correlation matrix) for the Gaussian copula.")
    }
    # Check that the correlation matrix is valid
    if (!is.matrix(copula_param) || nrow(copula_param) != cop_dim || ncol(copula_param) != cop_dim){
      stop("copula_param must be a square matrix with dimensions equal to the number of endpoints.")
    }
    # Ensure the correlation matrix is positive definite
    if (!all(eigen(copula_param)$values > 0)){
      stop("The provided correlation matrix is not positive definite.")
    }
    # Create the Gaussian copula
    cop <- normalCopula(param = P2p(copula_param), dim = cop_dim, dispstr = "un")
  } else {
    stop("Unsupported copula type. Please choose from 'Clayton', 'Frank', 'Gumbel', or 'Gaussian'.")
  }
  
  # Generate samples from the copula
  sample_copula <- rCopula(N.Super, cop)
  
  # Initialize lists to store data
  data_list <- vector("list", N.endpoints)
  survival_indices <- c()
  ordinal_indices <- c()
  binary_indices <- c()
  continuous_indices <- c()
  count_indices <- c()
  
  for (i in seq_along(endpoints)){
    endpoint <- endpoints[[i]]
    U <- sample_copula[, i] # Uniform(0,1) sample for this endpoint
    
    if (endpoint$type == "survival"){
      survival_indices <- c(survival_indices, i)
      dist <- endpoint$dist
      params <- endpoint$params
      # Generate survival times using the inverse CDF (quantile function)
      if (dist == "Exponential"){
        rate <- params$lambda
        T <- qexp(U, rate = rate)
      } else if (dist == "Weibull"){
        shape <- params$shape
        scale <- params$scale
        T <- qweibull(U, shape = shape, scale = scale)
      } else {
        stop(paste("Unsupported distribution for:", dist))
      }
      data_list[[i]] <- list(T = T)
    } else if (endpoint$type == "ordinal"){
      ordinal_indices <- c(ordinal_indices, i)
      prob <- endpoint$prob
      cumprob <- cumsum(prob)
      categories <- findInterval(U, cumprob) + 1 # Categories start from 1
      data_list[[i]] <- list(Categories = categories)
    } else if (endpoint$type == "binary"){
      binary_indices <- c(binary_indices, i)
      P <- endpoint$prob
      Y <- as.numeric(U <= P)
      data_list[[i]] <- list(Y = Y)
    } else if (endpoint$type == "continuous"){
      continuous_indices <- c(continuous_indices, i)
      mu <- endpoint$params$mu
      sigma <- endpoint$params$sigma
      Y <- qnorm(U, mean = mu, sd = sigma)
      data_list[[i]] <- list(Y = Y)
    } else if (endpoint$type == "count"){
      count_indices <- c(count_indices, i)
      lambda <- endpoint$params$lambda
      Y <- qpois(U, lambda = lambda)
      data_list[[i]] <- list(Y = Y)
    } else {
      stop(paste("Unsupported endpoint type:", endpoint$type))
    }
  }
  
  # Generate censoring times
  C <- rep(Fellow_up.Time, N.Super)
  
  # Initialize data frame
  data <- data.frame(Censoring_Time = C)
  
  # Process survival endpoints
  for (i in survival_indices){
    T <- data_list[[i]]$T
    # Apply censoring
    Y <- pmin(T, C)
    Delta <- as.numeric(T <= C)
    
    data[[paste0("T_", i)]] <- T
    data[[paste0("Y_", i)]] <- Y
    data[[paste0("delta_", i)]] <- Delta
  }
  
  # Process ordinal endpoints
  for (i in ordinal_indices){
    categories <- data_list[[i]]$Categories
    endpoint <- endpoints[[i]]
    data[[paste0("Ordinal_", i)]] <- factor(categories, levels = 1:length(endpoint$prob))
  }
  
  # Process binary endpoints
  for (i in binary_indices){
    Y <- data_list[[i]]$Y
    data[[paste0("Binary_", i)]] <- Y
  }
  
  # Process continuous endpoints
  for (i in continuous_indices){
    Y <- data_list[[i]]$Y
    data[[paste0("Continuous_", i)]] <- Y
  }
  
  # Process count endpoints
  for (i in count_indices){
    Y <- data_list[[i]]$Y
    data[[paste0("Count_", i)]] <- Y
  }
  
  return(data)
}
