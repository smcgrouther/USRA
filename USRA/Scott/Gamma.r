# General Gamma Distortion

gamma_integral <- function(a, b, c) {
    if (c == Inf) {
        return (beta(a,b))
        
    } else {
        integrand <- function(t) {
            t^(a - 1) * (1 - t)^(b - 1) * exp(-t / c)
        }
    }
    
    result <- integrate(integrand, lower = 0, upper = 1)
    
    return(result$value)
}


gamma_mean <- function(a, b, c, rho) {
    #' Function to calculate mean of gamma distorted distribution
    #'
    #' @param a - Parameter 'a' for gamma distortion.
    #' @param b - Parameter 'b' for gamma distortion.
    #' @param c - Parameter 'c' for gamma distortion.
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @returns Gamma_Mean - The mean of Z in the Vasicek Model.
    #' @use - Gamma_Mean <- gamma_mean(a, b, c, rho)
    #'

    gamma <- gamma_integral(a, b, c)
    # Define the integrand function 
    # z * 1/B(a,b) * PHI(z)^(a-1) * PHI(-z)^(b-1) * phi(z) ~ McGrouther ***
    integrand <- function(z, a, b, c, rho) {
    part1 <- z
    part2 <- 1 / gamma
    part3 <- pnorm(z)^(a - 1) * pnorm(-z)^(b - 1)
    part4 <- exp(-pnorm(z)/c)
    part5 <- dnorm(z)
    
    part1 * part2 * part3 * part4 * part5
        
    }
  
    # Integral of Above
    result <- integrate(integrand, lower = -37.01, upper = 37.01, a = a, b = b, c= c, rho = rho)

    # Obtaining The Values
    Gamma_Mean <- result$value
  
    return (Gamma_Mean)
    
}

gamma_variance <- function(a, b, c, rho) {
    #' Function to calculate variance of gamma distorted distribution
    #'
    #' @param a - Parameter 'a' for gamma distortion.
    #' @param b - Parameter 'b' for gamma distortion.
    #' @param c - Parameter 'c' for gamma distortion.
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @returns Gamma_Var - The varaince of Z in the Vasicek Model.
    #' @use - Gamma_Var <- gamma_variance(a, b, c, rho)
    #'

    # Calculate mean
    mean <- gamma_mean(a, b, c, rho)
  
    # Define the integrand function 
    # (z-mu)^2 * 1/B(a,b) * PHI(z)^(a-1) * PHI(-z)^(b-1) * phi(z) ~ McGrouther ***
    gamma <- gamma_integral(a, b, c)
    integrand <- function(z, a, b, c, rho) {
        part1 <- (z - mean)^2
        part2 <- 1 / gamma
        part3 <- pnorm(z)^(a - 1) * pnorm(-z)^(b - 1)
        part4 <- exp(-pnorm(z)/c)
        part5 <- dnorm(z)
    
        part1 * part2 * part3 * part4 * part5
  
    }
  
    # Integral of Above
    result <- integrate(integrand, lower = -37.01, upper = 37.01, a = a, b = b, c = c, rho = rho)

    # Obtaining The Values
    Gamma_Var <- result$value
  
    return (Gamma_Var)
    
}

gamma_default_prob <- function(a, b, c, PD, rho) {
    #' Function to calculate the distorted default probability using gamma 
    #'                       distortion. P(Y_n<= K_n)
    #'
    #' @param a - Parameter 'a' for gamma distortion.
    #' @param b - Parameter 'b' for gamma distortion.
    #' @param c - Parameter 'c' for gamma distortion.
    #' @param PD - Undistorted probability of default (Numeric value between 0 and 1).
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @returns Gamma_DDP - P(Y_n<=K_n) for Beta Distortion on Z ~ N(0,1)
    #' @use - Gamma_DDP <- gamma_default_prob(a, b, c, PD, rho)
    #'

    threshold <- qnorm(PD)  # Threshold
    
    # Define the integrand function 
    # PHI((K_n-sqrt(rho)*z) / sqrt(1-rho))* 1/B(a,b)* PHI(z)^(a-1)* PHI(-z)^(b-1)*phi(z) 
    # McGrouther ***

    gamma <- gamma_integral(a, b, c)
    integrand <- function(z, a, b, c, PD, rho) {
        part1 <- pnorm((threshold - sqrt(rho) * z) / sqrt(1 - rho))
        part2 <- 1 / gamma
        part3 <- pnorm(z)^(a - 1) * pnorm(-z)^(b - 1)
        part4 <- exp(-pnorm(z)/c)
        part5 <- dnorm(z)
    
        part1 * part2 * part3 * part4 * part5
    
    }
    
    # Integral of Above
    result <- integrate(integrand, lower = -37.01, upper = 37.01, a = a, b = b, c = c, PD = PD, rho = rho)

    # Obtaining The Values
    Gamma_DDP <- result$value
  
    return (Gamma_DDP)
    
}


gamma_default_corr <- function(a, b, c, PD, rho) {
    #' Function to calculate the distorted default correlation using gamma 
    #'                     distortion. Corr(I_D_n, I_D_m)
    #'
    #' @param a - Parameter 'a' for gamma distortion.
    #' @param b - Parameter 'b' for gamma distortion.
    #' @param c - Parameter 'c' for gamma distortion.
    #' @param PD - Undistorted probability of default (Numeric value between 0 and 1).
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @returns Gamma_DDC - Corr(I_D_n, I_D_m) (Numeric value between 0 and 1).
    #' @use - Gamma_DDC <- gamma_default_corr(a, b, c, PD, rho)
    #'
    
    threshold <- qnorm(PD)  # Threshold

    # Define the integrand function 
    # PHI((K_n-sqrt(rho)*z) / sqrt(1-rho))^2*1/B(a,b)*PHI(z)^(a-1)*PHI(-z)^(b-1)*phi(z) 
    # ~ McGrouther ***

    gamma <- gamma_integral(a, b, c)
    integrand <- function(z, a, b, c, PD, rho) {
        part1 <- (pnorm((threshold - sqrt(rho) * z) / sqrt(1 - rho)))^2
        part2 <- 1 / gamma
        part3 <- pnorm(z)^(a - 1) * pnorm(-z)^(b - 1)
        part4 <- exp(-pnorm(z)/c)
        part5 <- dnorm(z)
    
        value <- part1 * part2 * part3 * part4 * part5
    
        return(value)
      
    }

    # Integral of Above
    result1 <- integrate(integrand, lower = -37.01, upper = 37.01, a = a, b = b, c = c, PD = PD, rho = rho)
    
    # P(Y_n <= K_n)
    result2 <- gamma_default_prob(a, b, c, PD, rho)

    # Obtaining the Values
    value1 <- result1$value
    value2 <- result2

    # Formula
    gamma_DDC <- (value1 - value2^2) / (value2 * (1 - value2))
  
    return (gamma_DDC)
    
}


find_rho_gamma <- function(target_DDC, a, b, c, PD, tol = 1e-10) {
    #' Function to find the latent correlation rho given the target Gamma Distorted 
    #'                      Default Correlation (target_DDC)
    #' 
    #' @param target_DDC - target Gamma DDC (Numeric between 0 and 1).
    #' @param a - Parameter 'a' for gamma distortion.
    #' @param b - Parameter 'b' for gamma distortion.
    #' @param c - Parameter 'c' for gamma distortion.
    #' @param PD - the probability of default (Numeric value between 0 and 1).
    #' @param tol - tolerance for the root-finding algorithm (default is 1e-10).
    #' @return rho - the latent correlation between obligors (Numeric between 0 and 1).
    #' @use - rho <- find_rho_gamma(target_DDC, a, b, PD)
    #' 
  
    # Define the difference function
    difference_function <- function(rho) {
  
        calculated_DDC <- gamma_default_corr(a, b, c, PD, rho)
    
    return (calculated_DDC - target_DDC)
        
    }
  
    # Use uniroot to find the root of the difference function
    result <- uniroot(difference_function, lower = 0, upper = 1, tol = tol)
    rho <- result$root
  
    return (rho)
    
}