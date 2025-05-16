# Load necessary libraries
library(Matrix)

CALC_Kernal.Matrix <- function(Group.Treat, Group.Control, endpoints){
  M <- nrow(Group.Treat)
  N <- nrow(Group.Control)
  num_endpoints <- length(endpoints)
  
  # Initialize lists to store W_ij, L_ij, Omega_ij matrices for each endpoint
  Endpoint_Matrices <- vector("list", num_endpoints)
  
  # Initialize Omega_Kernal to ones for the starting point
  Omega_Kernal <- matrix(1, nrow = M, ncol = N)
  
  # Loop over endpoints according to their hierarchy
  for (i in seq_along(endpoints)){
    endpoint <- endpoints[[i]]
    endpoint_type <- endpoint$type
    endpoint_name <- paste0(endpoint$type, "_", i)
    
    # Initialize W_ij, L_ij, Omega_ij matrices
    W_ij <- Matrix(0, nrow = M, ncol = N, sparse = TRUE)
    L_ij <- Matrix(0, nrow = M, ncol = N, sparse = TRUE)
    Omega_ij <- Matrix(0, nrow = M, ncol = N, sparse = TRUE)
    
    if (endpoint_type == "survival"){
      # Extract survival times and event indicators
      Y_Treat <- Group.Treat[[paste0("Y_", i)]]
      Delta_Treat <- Group.Treat[[paste0("delta_", i)]]
      Y_Control <- Group.Control[[paste0("Y_", i)]]
      Delta_Control <- Group.Control[[paste0("delta_", i)]]
      # Calculate W_ij and L_ij matrices
      W_ij <- (outer(Y_Treat, Y_Control, ">") & matrix(Delta_Control, nrow = M, ncol = N, byrow = TRUE)) * 1
      L_ij <- (outer(Y_Treat, Y_Control, "<") & matrix(Delta_Treat, nrow = M, ncol = N, byrow = FALSE)) * 1
      Omega_ij <- (1 - W_ij - L_ij)
    } 
    else if (endpoint_type == "ordinal"){
      # Extract ordinal values
      Y_Treat <- as.integer(Group.Treat[[paste0("Ordinal_", i)]])
      Y_Control <- as.integer(Group.Control[[paste0("Ordinal_", i)]])
      # For ordinal variables, lower values are better
      W_ij <- (outer(Y_Treat, Y_Control, ">")) * 1
      L_ij <- (outer(Y_Treat, Y_Control, "<")) * 1
      Omega_ij <- (outer(Y_Treat, Y_Control, "==")) * 1
    } 
    else if (endpoint_type == "binary"){
      # Extract binary values
      Y_Treat <- Group.Treat[[paste0("Binary_", i)]]
      Y_Control <- Group.Control[[paste0("Binary_", i)]]
      # Assuming Y = 1 is success, Y = 0 is event
      W_ij <- (outer(Y_Treat, Y_Control, ">")) * 1
      L_ij <- (outer(Y_Treat, Y_Control, "<")) * 1
      Omega_ij <- (outer(Y_Treat, Y_Control, "==")) * 1
    } 
    else if (endpoint_type == "continuous"){
      # Extract continuous values
      Y_Treat <- Group.Treat[[paste0("Continuous_", i)]]
      Y_Control <- Group.Control[[paste0("Continuous_", i)]]
      # Use the threshold provided in the endpoint specification
      threshold <- endpoint$threshold
      # W_ij is 1 when (Y_Treat - Y_Control) > threshold
      W_ij <- outer(Y_Treat, Y_Control, FUN = function(x, y) { (x - y) > threshold }) * 1
      # L_ij is 1 when (Y_Control - Y_Treat) > threshold
      L_ij <- outer(Y_Treat, Y_Control, FUN = function(x, y) { (y - x) > threshold }) * 1
      Omega_ij <- 1 - W_ij - L_ij
    } 
    else if (endpoint_type == "count"){
      # Extract count values
      Y_Treat <- Group.Treat[[paste0("Count_", i)]]
      Y_Control <- Group.Control[[paste0("Count_", i)]]
      # Assuming lower counts are better (e.g., fewer adverse events)
      W_ij <- (outer(Y_Treat, Y_Control, "<")) * 1
      L_ij <- (outer(Y_Treat, Y_Control, ">")) * 1
      Omega_ij <- (outer(Y_Treat, Y_Control, "==")) * 1
    } 
    else {
      stop(paste("Unsupported endpoint type:", endpoint_type))
    }
    # Multiply by previous Omega_Kernal to account for hierarchy
    W_ij <- Omega_Kernal * W_ij
    L_ij <- Omega_Kernal * L_ij
    Omega_ij <- Omega_Kernal * Omega_ij
    
    # Store the matrices for endpoint i
    Endpoint_Matrices[[i]] <- list( W_ij = W_ij, L_ij = L_ij, Omega_ij = Omega_ij )
    
    # Update Omega_Kernal for next endpoint
    Omega_Kernal <- Omega_ij
  }
  
  # Sum the W_ij and L_ij matrices across endpoints to get Win_Kernal and Loss_Kernal
  Win_Kernal <- Reduce("+", lapply(Endpoint_Matrices, function(x) x$W_ij))
  Loss_Kernal <- Reduce("+", lapply(Endpoint_Matrices, function(x) x$L_ij))
  
  # Calculate win and loss probabilities
  tau_w <- mean(as.vector(Win_Kernal > 0))
  tau_l <- mean(as.vector(Loss_Kernal > 0))
  
  # You can also extract individual tau_w and tau_l for each endpoint if needed
  tau_w_list <- sapply(Endpoint_Matrices, function(x) mean(as.vector(x$W_ij)))
  tau_l_list <- sapply(Endpoint_Matrices, function(x) mean(as.vector(x$L_ij)))
  
  return(list(
    Win_Kernal = Win_Kernal,
    Loss_Kernal = Loss_Kernal,
    tau_w = tau_w,
    tau_l = tau_l,
    tau_w_list = tau_w_list,
    tau_l_list = tau_l_list
  ))
}

CALC_Xi <- function(Win_Kernal, Loss_Kernal){
  M <- nrow(Win_Kernal); N <- ncol(Win_Kernal)
  tau_w <- mean(Win_Kernal); tau_l <- mean(Loss_Kernal)
  
  xi.ww10 <- (1 / (M * N * (N - 1))) * (sum(rowSums(Win_Kernal) * rowSums(Win_Kernal)) - sum(Win_Kernal * Win_Kernal)) - tau_w^2
  xi.wl10 <- (1 / (M * N * (N - 1))) * (sum(rowSums(Win_Kernal) * rowSums(Loss_Kernal)) - sum(Win_Kernal * Loss_Kernal)) - tau_w * tau_l
  xi.ll10 <- (1 / (M * N * (N - 1))) * (sum(rowSums(Loss_Kernal) * rowSums(Loss_Kernal)) - sum(Loss_Kernal * Loss_Kernal)) - tau_l^2
  
  xi.ww01 <- (1 / (M * N * (M - 1))) * (sum(colSums(Win_Kernal) * colSums(Win_Kernal)) - sum(Win_Kernal * Win_Kernal)) - tau_w^2
  xi.wl01 <- (1 / (M * N * (M - 1))) * (sum(colSums(Win_Kernal) * colSums(Loss_Kernal)) - sum(Win_Kernal * Loss_Kernal)) - tau_w * tau_l
  xi.ll01 <- (1 / (M * N * (M - 1))) * (sum(colSums(Loss_Kernal) * colSums(Loss_Kernal)) - sum(Loss_Kernal * Loss_Kernal)) - tau_l^2
  
  xi.ww11 <- 1 / (M * N) * (sum(Win_Kernal * Win_Kernal)) - tau_w^2
  xi.wl11 <- 1 / (M * N) * (sum(Win_Kernal * Loss_Kernal))  - tau_w * tau_l
  xi.ll11 <- 1 / (M * N) * (sum(Loss_Kernal * Loss_Kernal)) - tau_l^2
  
  return(data.frame(xi.ww10 = xi.ww10, xi.wl10 = xi.wl10, xi.ll10 = xi.ll10, 
                    xi.ww01 = xi.ww01, xi.wl01 = xi.wl01, xi.ll01 = xi.ll01,
                    xi.ww11 = xi.ww11, xi.wl11 = xi.wl11, xi.ll11 = xi.ll11))
}


#### Variance of a series of derived Win Statistics: (NB, WR, WO and DOOR)
Var_NB<- function(m, n, Xi){
  Var.NB <- (n - 1) / (m * n) * (Xi$xi.ww10 + Xi$xi.ll10 - 2 * Xi$xi.wl10) + (m - 1) / (m * n) * (Xi$xi.ww01 + Xi$xi.ll01 - 2 * Xi$xi.wl01) + 1 / (m * n) * (Xi$xi.ww11 + Xi$xi.ll11 - 2 * Xi$xi.wl11)
  return(Var.NB)
}

Var_logWR <- function(m, n, Xi, tau_w, tau_l){
  Var.logWR <- ((n - 1) / (m * n) * Xi$xi.ww10 + (m - 1) / (m * n) * Xi$xi.ww01 + 1 / (m * n) * Xi$xi.ww11) / tau_w^2 +
    ((n - 1) / (m * n) * Xi$xi.ll10 + (m - 1) / (m * n) * Xi$xi.ll01 + 1 / (m * n) * Xi$xi.ll11) / tau_l^2  -
    2 *  ((n - 1) / (m * n) * Xi$xi.wl10 + (m - 1) / (m * n) * Xi$xi.wl01 + 1 / (m * n) * Xi$xi.wl11)/ (tau_w * tau_l)
  return(Var.logWR)
}

Var_logWO <- function(m, n, Xi, tau_w, tau_l){
  Var.NB <- Var_NB(m, n, Xi)
  delta.tau <- tau_w - tau_l
  Var.logWO <- 4 * Var.NB / (delta.tau^2 -1)^2
  return(Var.logWO)
}

Var_DOOR <- function(m, n, Xi){
  Var.NB <- Var_NB(m, n, Xi)
  return(1/4 * Var.NB)
}

# Sample size calculation function
SampleSize_Calc <- function(tau_w.HA, tau_l.HA, Xi.HA, Xi.H0, tau_w.H0, tau_l.H0, 
                            alpha, beta, Sample.rho, Metric){
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(1 - beta)
  # Calculate delta_tau and set up variance functions based on the Metric
  if(Metric == "NB"){
    delta_tau <- tau_w.HA - tau_l.HA
    f <- function(m){    # Define the function to find the root
      n <- m * Sample.rho
      Var_H0 <- Var_NB(m = m, n = n, Xi = Xi.H0)
      Var_HA <- Var_NB(m = m, n = n, Xi = Xi.HA)
      lhs <- - z_beta * sqrt(Var_HA) # z_{1-\beta} = - z_\beta
      rhs <- z_alpha * sqrt(Var_H0) - delta_tau
      return(lhs - rhs)
    }
  } else if(Metric == "WR"){
    theta_tau <- log(tau_w.HA / tau_l.HA)
    f <- function(m){
      n <- m * Sample.rho
      Var_H0 <- Var_logWR(m = m, n = n, Xi = Xi.H0, tau_w = tau_w.H0, tau_l = tau_l.H0)
      Var_HA <- Var_logWR(m = m, n = n, Xi = Xi.HA, tau_w = tau_w.HA, tau_l = tau_l.HA)
      lhs <- - z_beta * sqrt(Var_HA) 
      rhs <- z_alpha * sqrt(Var_H0) - theta_tau
      return(lhs - rhs)
    }
  } else if(Metric == "WO"){
    gamma_tau <- log(0.5*(1 + tau_w.HA - tau_l.HA)) - log(0.5*(1 - tau_w.HA + tau_l.HA))
    # Define the function to find the root
    f <- function(m){
      n <- m * Sample.rho
      Var_H0 <- Var_logWO(m = m, n = n, Xi = Xi.H0, tau_w = tau_w.H0, tau_l = tau_l.H0)
      Var_HA <- Var_logWO(m = m, n = n, Xi = Xi.HA, tau_w = tau_w.HA, tau_l = tau_l.HA)
      lhs <- -z_beta * sqrt(Var_HA)
      rhs <- z_alpha * sqrt(Var_H0) - gamma_tau
      return(lhs - rhs)
    }
  } else if(Metric == "DOOR"){
    lambda_tau <- 0.5*(1 + tau_w.HA - tau_l.HA)
    # Define the function to find the root
    f <- function(m){
      n <- m * Sample.rho
      Var_H0 <- Var_DOOR(m = m, n = n, Xi = Xi.H0)
      Var_HA <- Var_DOOR(m = m, n = n, Xi = Xi.HA)
      lhs <- - z_beta * sqrt(Var_HA)
      rhs <- z_alpha * sqrt(Var_H0) - lambda_tau + 0.5
      return(lhs - rhs)
    }
  } else {
    stop("Invalid Metric. Choose one of 'NB', 'WR', 'WO', 'DOOR'.")
  }
  
  # Use uniroot to solve for m
  lower_bound <- 10
  upper_bound <- 1e6
  
  f_lower <- f(lower_bound)
  f_upper <- f(upper_bound)
  
  if (f_lower * f_upper > 0){
    stop("Root not found in the specified interval. Please adjust the bounds.")
  }
  
  result <- uniroot(f, lower = lower_bound, upper = upper_bound)
  
  m.sample<- ceiling(result$root)
  n.sample <- ceiling(m.sample * Sample.rho)
  
  return(list(
    m.sample = m.sample,
    n.sample = n.sample
  ))
}

# Power calculation function (Given Sample size) return theortical power level
# m is the sample size for treatment group and sample size in ctrl n can be calculated as n = m * Sample.rho
Power_Calc <- function(tau_w.HA, tau_l.HA, Xi.HA, Xi.H0, tau_w.H0, tau_l.H0, 
                       alpha, m, Sample.rho, Metric){
  z_alpha <- qnorm(1 - alpha)
  n <- m * Sample.rho
  # Calculate delta_tau and set up variance functions based on the Metric
  if(Metric == "NB"){
    delta_tau <- tau_w.HA - tau_l.HA
    Var_H0 <- Var_NB(m = m, n = n, Xi = Xi.H0)
    Var_HA <- Var_NB(m = m, n = n, Xi = Xi.HA)
    power <- pnorm(q = (z_alpha * sqrt(Var_H0) - delta_tau) / sqrt(Var_HA), lower.tail = FALSE )
  }
  else if(Metric == "WR"){
    theta_tau <- log(tau_w.HA / tau_l.HA)
    Var_H0 <- Var_logWR(m = m, n = n, Xi = Xi.H0, tau_w = tau_w.H0, tau_l = tau_l.H0)
    Var_HA <- Var_logWR(m = m, n = n, Xi = Xi.HA, tau_w = tau_w.HA, tau_l = tau_l.HA)
    power <- pnorm(q = (z_alpha * sqrt(Var_H0) - theta_tau) / sqrt(Var_HA), lower.tail = FALSE )
  } 
  else if(Metric == "WO"){
    gamma_tau <- log(0.5*(1 + tau_w.HA - tau_l.HA)) - log(0.5*(1 - tau_w.HA + tau_l.HA))
    Var_H0 <- Var_logWO(m = m, n = n, Xi = Xi.H0, tau_w = tau_w.H0, tau_l = tau_l.H0)
    Var_HA <- Var_logWO(m = m, n = n, Xi = Xi.HA, tau_w = tau_w.HA, tau_l = tau_l.HA)
    power <- pnorm(q = (z_alpha * sqrt(Var_H0) - gamma_tau) / sqrt(Var_HA), lower.tail = FALSE )
  } 
  else if(Metric == "DOOR"){
    lambda_tau <- 0.5*(1 + tau_w.HA - tau_l.HA)
    Var_H0 <- Var_DOOR(m = m, n = n, Xi = Xi.H0)
    Var_HA <- Var_DOOR(m = m, n = n, Xi = Xi.HA)
    power <- pnorm(q = (z_alpha * sqrt(Var_H0) - - lambda_tau + 0.5) / sqrt(Var_HA), lower.tail = FALSE )
  } 
  else {
    stop("Invalid Metric. Choose one of 'NB', 'WR', 'WO', 'DOOR'.")
  }
  return(power)
}

