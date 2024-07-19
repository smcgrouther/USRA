# Building the Loss Distribution, Using the Vasicek Model

single_threshold <- function(n, rho) {
    #' Function to generate random credit scores based on the Vasicek Model.
    #' 
    #' @param n - Amount of obligors (Positive integer).
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @return Y - Credit scores of n obligors (Vector of size n).
    #' 
    #' @use Y <- single_threshold(n, rho)
    
    U <- runif(1)
    Z <- qnorm(U)
  
    U2 <- runif(n)
    epsilon <- qnorm(U2)
  
    Y <- sqrt(rho) * Z + sqrt(1 - rho) * epsilon
  
    return(Y)
    
}

distorted_LC <- function(rho, variance) {
    #' Function to calculate the distorted latent correlation given the variance of Z.
    #' 
    #' @param rho - Latent correlation (Numeric value between 0 and 1).
    #' @param variance - Variance of Z in the Vasicek Model Model (Positive real number).
    #' @return rho_d - Distorted latent correlation (Numeric value between 0 and 1).
    #'
    #' @use rho_d <- distorted_LC(rho, variance)
    
    rho_d <- rho * variance / (rho * variance + 1 - rho)

    return(rho_d)
    
}

default_amount <- function(n, rho, PD, mu, sigma, a, b) {
    #' Function to generate which obligors defaulted using the Vasicek Model.
    #' 
    #' @param n - Amount of obligors (Positive integer).
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @param PD - Undistorted probability of default (Numeric value between 0 and 1).
    #' @param mu - Mean of gaussian distortion (Numeric value).
    #' @param sigma - Standard deviation of gaussian distortion (Numeric value > 0).
    #' @param a - Parameter 'a' for beta distortion (Optional, defaults to 1).
    #' @param b - Parameter 'b' for beta distortion (Optional, defaults to 1).
    #' @return default_list - Default indication under each distortion type (nx3 matrix).
    #' @use default_list <- default_amount(n, rho, PD, mu, sigma, a, b)
    #' Notes: Col 1 => Simple, Col 2 => Gaussian, Col 3 => Beta
    #' 

    # Thresholds
    threshold_simple <- qnorm(PD)
    threshold_gaussian <- qnorm(gaussian_default_prob(mu, sigma, PD, rho))
    threshold_beta <- qnorm(beta_default_prob(a, b, PD, rho))

    # Credit Scores
    Y <- single_threshold(n, rho)

    # Definitions
    default_simple <- numeric(n)
    default_gaussian <- numeric(n)
    default_beta <- numeric(n)

    # Determining if Default Occours Under Each Model
    for (i in 1:n) {
        default_simple[i] <- ifelse(Y[i] <= threshold_simple, 1, 0)
        default_gaussian[i] <- ifelse(Y[i] <= threshold_gaussian, 1, 0)
        default_beta[i] <- ifelse(Y[i] <= threshold_beta, 1, 0)
    }

    # Combining All Into a Matrix
    default_list <- cbind(default_simple, default_gaussian, default_beta)
  
    return(default_list)
    
}

loss_generation <- function(EAD, RR, n, rho, PD, mu, sigma, a, b) {
    #'           Function to generate the loss of each obligor 
    #'          under each distortion type using the Vasicek Model.
    #' 
    #' @param EAD - Exposure at Default (Numeric value).
    #' @param RR - Recovery Rate (Numeric value between 0 and 1).
    #' @param n - Amount of obligors (Positive integer).
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @param PD - Undistorted probability of default (Numeric value between 0 and 1).
    #' @param mu - Mean of gaussian distortion (Numeric value).
    #' @param sigma - Standard deviation of gaussian distortion (Numeric value > 0).
    #' @param a - Parameter 'a' for beta distortion.
    #' @param b - Parameter 'b' for beta distortion.
    #' @return loss - Loss per customer under each distortion type (nx3 matrix).  
    #' @use loss <- loss_generation(EAD, RR, n, rho, PD, mu, sigma, a, b)
    #' Notes: Col 1 => Simple, Col 2 => Gaussian, Col 3 => Beta
    #'

    # Function to generate which obligors defaulted using the Vasicek Model.
    default_list <- default_amount(n, rho, PD, mu, sigma, a, b)

    # Loss Formula
    loss <- EAD * (1 - RR) * default_list   

    return(loss)
    
}

total_loss <- function(loss) {
    #' Function to generate the total loss under each distortion using the Vasicek Model.
    #' 
    #' @param loss - Loss per customer under each distortion type (nx3 matrix).
    #' @return total - Total loss under each distortion type (1x3 vector).
    #' @use total <- total_loss(loss)
    #' 
    
    total <- colSums(loss)
    
    return(total)
    
}

monte_carlo_loss <- function(MC, EAD, RR, n, rho, PD, mu, sigma, a, b) {
    #'     Monte Carlo simulation to estimate the loss per customer under each distortion 
    #'              type using the Single Factor Threshold Model (Vasicek Model)
    #'
    #' @param MC - Amount of simulations (Positive integer).
    #' @param EAD - Exposure at Default (Numeric value).
    #' @param RR - Recovery Rate (Numeric value between 0 and 1).
    #' @param n - Amount of obligors (Positive integer).
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @param PD - Undistorted probability of default (Numeric value between 0 and 1).
    #' @param mu - Mean of gaussian distortion (Numeric value).
    #' @param sigma - Standard deviation of gaussian distortion (Numeric value > 0).
    #' @param a - Parameter 'a' for beta distortion.
    #' @param b - Parameter 'b' for beta distortion.
    #' @return mc_loss - Matrix of total losses under each distortion type (MCx3).
    #' @use mc_loss <- monte_carlo_loss(MC, EAD, RR, n, rho, PD, mu, sigma, a, b)
    #' Notes:
    #'

    # Initalize
    mc_loss <- matrix(0, nrow = MC, ncol = 3)

    # Monte-Carlo Simulation
    for (i in 1:MC) {
        loss <- loss_generation(EAD, RR, n, rho, PD, mu, sigma, a, b)
        total <- total_loss(loss)
        mc_loss[i, ] <- total
    }
    
    return(mc_loss)
    
}

average_loss <- function(mc_loss) {
    #'       Function to calculate the average loss under each distortion type  
    #'            Given a Monte Carlo simulation using the Vasicek Model.
    #'
    #' @param mc_loss - Monte carlo simulation of total losses under each distortion type (MCx3 vector).
    #' @return avg - Average loss given a Monte-Carlo simulation (1x3 vector).
    #' @use avg <- average_loss(mc_loss)
    #'
    
    avg <- colMeans(mc_loss)
    
    return(avg)
    
}

value_at_risk <- function(mc_loss, alpha) {
    #' Function to calculate the value at risk at the alpha rate through a Monte Carlo 
    #'                      simulation using the Vasicek Model.
    #' 
    #' @param mc_loss - Monte carlo simulation of total losses under each distortion type (MCx3 vector).
    #' @param alpha - Alpha rate (Numeric value between 0 and 1).
    #' @return VaR - Value at risk at the alpha rate for each distortion type (1x3 vector).
    #' @use VaR <- value_at_risk(mc_loss, alpha)
    #' 
    
    sorted_mc_loss <- apply(mc_loss, 2, sort)
    size <- nrow(mc_loss)
    at_value <- ceiling(size * (1-alpha))
  
    VaR <- sorted_mc_loss[at_value, ]
  
    return(VaR)
    
}

expected_shortfall <- function(mc_loss, alpha) {
    #' Function to calculate the expected shortfall at the alpha rate through a Monte 
    #'                     Carlo simulation using the Vasicek Model.                      
    #' 
    #' @param mc_loss - Monte carlo simulation of total losses under each distortion type (MCx3 vector).
    #' @param alpha - Alpha rate (Numeric value between 0 and 1).
    #' @return Expected shortfall at the alpha rate for each distortion type (1x3 vector).
    #' @use ES <- expected_shortfall(mc_loss, alpha)
    #'
    
    sorted_mc_loss <- apply(mc_loss, 2, sort)
    size <- nrow(mc_loss)
    at_value <- ceiling(size * (1-alpha))
  
    ES <- colMeans(sorted_mc_loss[at_value:size, ,  drop = FALSE])
  
    return(ES)
    
}

confidence_intervals <- function(c, MC, EAD, RR, n, rho, PD, mu, sigma, a, b, alpha) {
    AVG <- matrix(0, nrow = c, ncol = 3)
    VAR <- array(0, dim = c(length(alpha), 3, c))
    ES <- array(0, dim = c(length(alpha), 3, c))

    for (i in 1:c) {
        mc_loss <- monte_carlo_loss(MC, EAD, RR, n, rho, PD, mu, sigma, a, b)
        AVG[i, ] <- average_loss(mc_loss)

        for (j in 1:length(alpha)) {
            VAR[j, , i] <- value_at_risk(mc_loss, alpha[j])
            ES[j, , i] <- expected_shortfall(mc_loss, alpha[j])
        }
    }

    avg_results <- colMeans(AVG)  # Calculate the average across all simulations

    VAR_sum <- apply(VAR, c(1, 2), sum)
    ES_sum <- apply(ES, c(1, 2), sum)

    return(list(ci_avg = avg_results, ci_VaR = VAR_sum, ci_ES = ES_sum))
}