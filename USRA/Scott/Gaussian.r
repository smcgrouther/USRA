# GAUSSIAN DISTORTION FUNCTIONS

gaussian_default_prob <- function(mu, sigma, PD, rho, bound = 37.01) {
    #' Function to generate the probability that an obligor defaulted under Gaussian 
    #'                         Distortion - P(Y_n <= K_n)
    #' @param mu - mean of Gaussian distortion (Numeric value).
    #' @param sigma - standard deviation of Gaussian distortion (Numeric value).
    #' @param PD - the probability of default (Numeric between 0 and 1).
    #' @param rho - the latent correlation between obligors (Numeric between 0 and 1).
    #' @return Gaussian_DDP - P(Y_n <= K_n) for Gaussian Distortion on Z ~ N(0,1)
    #' Use Gaussian_DDP <- gaussian_default_prob(mu, sigma, PD, rho)
  
    threshold <- qnorm(PD)  # Threshold

    # Define the integrand function 
    # PHI((K_n-sqrt(rho)*z)/sqrt(1-rho)) * 1/sigma * phi((z-mu)/sigma) ~ McGrouther ***
    integrand <- function(z, mu, sigma, PD, rho) {
        part1 <- pnorm((threshold - sqrt(rho) * z) / sqrt(1 - rho)) 
        part2 <- dnorm((z - mu) / sigma)                            
        part3 <- (1 / sigma)                                      
        part1 * part2 * part3
    }

    # Integral of Above
    result <- integrate(integrand, lower = -bound, upper = bound, mu = mu, sigma = sigma, PD = PD, rho = rho)

    # Obtaining The Values
    Gaussian_DDP <- result$value
  
    return (Gaussian_DDP)
  
}

# Gaussian Distorted Default Correlation - Corr(I_D_{n},I_D_{m})

gaussian_default_corr <- function(mu, sigma, PD, rho, bound = 37.01){
    #' Function to generate the default correlation.  Correlation that two similar 
    #'                   obligors defaulted using Gaussian Distortion.
    #' @param mu - mean of Gaussian distortion (Numeric value).
    #' @param sigma - standard deviation of Gaussian distortion (Numeric value).
    #' @param PD - the probability of default (Numeric between 0 and 1).
    #' @param rho - the latent correlation between obligors (Numeric between 0 and 1).
    #' @return Gaussian_DDC - corr(I_D{n}, I_D{m})
    #' Use Gaussian_DDC <- gaussian_default_corr(mu, sigma, PD, rho)
   
    threshold <- qnorm(PD)  # Threshold

    # Define the integrand function 
    # PHI((K_n-sqrt(rho)*z)/sqrt(1-rho))^2 * 1/sigma * phi((z-mu)/sigma) ~ McGrouther ***

    
    integrand <- function(z, mu, sigma, PD, rho) {
        part1 <- pnorm((threshold - sqrt(rho) * z) / sqrt(1 - rho))^2  
        part2 <- dnorm((z - mu) / sigma)                                
        part3 <- (1 / sigma)                                         
        part1 * part2 * part3
    }

    # Integral of Above

    result1 <- integrate(integrand, lower = -bound , upper = bound, mu = mu, sigma = sigma, PD = PD, rho = rho)

    
    # P(Y_n <= K_n)
    result2 <- gaussian_default_prob(mu, sigma, PD, rho, bound) 

    # Obtaining the Values
    value1 <- result1$value
    value2 <- result2

    # Formula
    Gaussian_DDC <- (value1 - value2^2)/((value2)*(1-value2))       
    
    return (Gaussian_DDC)
  
}

find_rho_gaussian <- function(target_DDC, mu, sigma, PD, tol = 1e-6) {
    #' Function to find the latent correlation rho given the target Gaussian Distorted 
    #' Default Correlation (Gaussian_DDC)
    #' @param target_DDC - target Gaussian DDC (Numeric between 0 and 1).
    #' @param mu - mean of Gaussian distortion (Numeric value).
    #' @param sigma - standard deviation of Gaussian distortion (Numeric value).
    #' @param PD - the probability of default (Numeric between 0 and 1).
    #' @param tol - tolerance for the root-finding algorithm (default is 1e-6).
    #' @return rho - the latent correlation between obligors (Numeric between 0 and 1).
    #' Use rho <- find_rho_gaussian(target_DDC, mu, sigma, PD)
  
    # Define the difference function
    difference_function <- function(rho) {
        calculated_DDC <- gaussian_default_corr(mu, sigma, PD, rho)
        return (calculated_DDC - target_DDC)
    }
  
    # Attempt to find root with initial bounds 0 to 1
    result <- try(uniroot(difference_function, lower = 0, upper = 1, tol = tol), silent = TRUE)

    # Check if uniroot failed or returned NA
    if (inherits(result, "try-error") || is.na(result$root)) {
        print("Root finding with initial bounds failed. Trying with reduced integration bounds.  May need to change manually within the code as it is hard fixed bounds if there is still an error.")
        
        # Define difference function with reduced integration bounds
        difference_function_reduced <- function(rho) {
            calculated_DDC <- gaussian_default_corr(mu, sigma, PD, rho, bound = 20)
            return (calculated_DDC - target_DDC)
        }
        
        # Try uniroot again with reduced integration bounds
        result <- uniroot(difference_function_reduced, lower = 0, upper = 1, tol = tol)
    }

    rho <- result$root
  
    return (rho)
}
