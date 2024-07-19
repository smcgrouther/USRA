# BETA DISTORTION FUNCTIONS

# Special Case (Beta Distortion)

beta_mean <- function(a, b, rho) {
    #' Function to calculate mean of beta/gamma distorted distribution
    #'
    #' @param a - Parameter 'a' for beta distortion.
    #' @param b - Parameter 'b' for beta distortion.
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @returns Beta_Mean - The mean of Z in the Vasicek Model.
    #' @use - Beta_Mean <- beta_mean(a, b, rho)
    #'

    # Define the integrand function 
    # z * 1/B(a,b) * PHI(z)^(a-1) * PHI(-z)^(b-1) * phi(z) ~ McGrouther ***
    integrand <- function(z, a, b, rho) {
    part1 <- z
    part2 <- 1 / beta(a, b)
    part3 <- pnorm(z)^(a - 1) * pnorm(-z)^(b - 1)
    part4 <- dnorm(z)
    
    part1 * part2 * part3 * part4
        
    }
  
    # Integral of Above
    result <- integrate(integrand, lower = -37.01, upper = 37.01, a = a, b = b, rho = rho)

    # Obtaining The Values
    Beta_Mean <- result$value
  
    return (Beta_Mean)
    
}

beta_variance <- function(a, b, rho) {
    #' Function to calculate variance of beta/gamma distorted distribution
    #'
    #' @param a - Parameter 'a' for beta distortion.
    #' @param b - Parameter 'b' for beta distortion.
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @returns Beta_Var - The varaince of Z in the Vasicek Model.
    #' @use - Beta_Var <- beta_variance(a, b, rho)
    #'

    # Calculate mean
    mean <- beta_mean(a, b, rho)
  
    # Define the integrand function 
    # (z-mu)^2 * 1/B(a,b) * PHI(z)^(a-1) * PHI(-z)^(b-1) * phi(z) ~ McGrouther ***
    integrand <- function(z, a, b, rho) {
        part1 <- (z - mean)^2
        part2 <- 1 / beta(a, b)
        part3 <- pnorm(z)^(a - 1) * pnorm(-z)^(b - 1)
        part4 <- dnorm(z)
    
        part1 * part2 * part3 * part4
  
    }
  
    # Integral of Above
    result <- integrate(integrand, lower = -37.01, upper = 37.01, a = a, b = b, rho = rho)

    # Obtaining The Values
    Beta_Var <- result$value
  
    return (Beta_Var)
    
}

beta_default_prob <- function(a, b, PD, rho) {
    #' Function to calculate the distorted default probability using beta/gamma 
    #'                       distortion. P(Y_n<= K_n)
    #'
    #' @param a - Parameter 'a' for beta distortion.
    #' @param b - Parameter 'b' for beta distortion.
    #' @param PD - Undistorted probability of default (Numeric value between 0 and 1).
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @returns Beta_DDP - P(Y_n<=K_n) for Beta Distortion on Z ~ N(0,1)
    #' @use - Beta_DDP <- beta_default_prob(a, b, PD, rho)
    #'

    threshold <- qnorm(PD)  # Threshold
    
    # Define the integrand function 
    # PHI((K_n-sqrt(rho)*z) / sqrt(1-rho))* 1/B(a,b)* PHI(z)^(a-1)* PHI(-z)^(b-1)*phi(z) 
    # McGrouther ***
    
    integrand <- function(z, a, b, PD, rho) {
        part1 <- pnorm((threshold - sqrt(rho) * z) / sqrt(1 - rho))
        part2 <- 1 / beta(a, b)
        part3 <- pnorm(z)^(a - 1) * pnorm(-z)^(b - 1)
        part4 <- dnorm(z)
    
        part1 * part2 * part3 * part4
    
    }
    
    # Integral of Above
    result <- integrate(integrand, lower = -37.01, upper = 37.01, a = a, b = b, PD = PD, rho = rho)

    # Obtaining The Values
    Beta_DDP <- result$value
  
    return (Beta_DDP)
    
}


beta_default_corr <- function(a, b, PD, rho) {
    #' Function to calculate the distorted default correlation using beta/gamma 
    #'                     distortion. Corr(I_D_n, I_D_m)
    #'
    #' @param a - Parameter 'a' for beta distortion.
    #' @param b - Parameter 'b' for beta distortion.
    #' @param PD - Undistorted probability of default (Numeric value between 0 and 1).
    #' @param rho - Latent correlation between obligors (Numeric value between 0 and 1).
    #' @returns Beta_DDC - Corr(I_D_n,I_D_m) (Numeric value between 0 and 1).
    #' @use - Beta_DDC <- beta_default_corr(a, b, PD, rho)
    #'
    
    threshold <- qnorm(PD)  # Threshold

    # Define the integrand function 
    # PHI((K_n-sqrt(rho)*z) / sqrt(1-rho))^2*1/B(a,b)*PHI(z)^(a-1)*PHI(-z)^(b-1)*phi(z) 
    # ~ McGrouther ***
    
    integrand <- function(z, a, b, PD, rho) {
        part1 <- pnorm((threshold - sqrt(rho) * z) / sqrt(1 - rho))^2
        part2 <- 1 / beta(a, b)
        part3 <- pnorm(z)^(a - 1) * pnorm(-z)^(b - 1)
        part4 <- dnorm(z)
    
        value <- part1 * part2 * part3 * part4
    
        return(value)
      
    }

    # Integral of Above
    result1 <- integrate(integrand, lower = -37.01, upper = 37.01, a = a, b = b, PD = PD, rho = rho)
    
    # P(Y_n <= K_n)
    result2 <- beta_default_prob(a, b, PD, rho)

    # Obtaining the Values
    value1 <- result1$value
    value2 <- result2

    # Formula
    Beta_DDC <- (value1 - value2^2) / (value2 * (1 - value2))
  
    return (Beta_DDC)
    
}


find_rho_beta <- function(target_DDC, a, b, PD, tol = 1e-10) {
    #' Function to find the latent correlation rho given the target Beta Distorted 
    #'                      Default Correlation (Beta_DDC)
    #' 
    #' @param target_DDC - target Beta DDC (Numeric between 0 and 1).
    #' @param a - Parameter 'a' for beta distortion.
    #' @param b - Parameter 'b' for beta distortion.
    #' @param PD - the probability of default (Numeric value between 0 and 1).
    #' @param tol - tolerance for the root-finding algorithm (default is 1e-10).
    #' @return rho - the latent correlation between obligors (Numeric between 0 and 1).
    #' @use - rho <- find_rho_beta(target_DDC, a, b, PD)
    #' 
  
    # Define the difference function
    difference_function <- function(rho) {
  
        calculated_DDC <- beta_default_corr(a, b, PD, rho)
    
    return (calculated_DDC - target_DDC)
        
    }
  
    # Use uniroot to find the root of the difference function
    result <- uniroot(difference_function, lower = 0, upper = 1, tol = tol)
    rho <- result$root
  
    return (rho)
    
}

# Proportional Hazards Distortion - Speacial Case

ph_mean <- function(a,rho){

    return (beta_mean(a,1,rho))

}

ph_variance <- function(a,rho){

    return(beta_variance(a,1,rho))

}

ph_default_prob <- function(a,PD,rho){

    return(beta_default_prob(a, 1, PD, rho))

}

ph_default_corr <- function(a,PD,rho){

    return(beta_default_corr(a, 1, PD, rho))

}

find_rho_ph <- function(target_DDC, a, PD, tol = 1e-10){

    return(find_rho_beta(target_DDC, a, 1, PD, tol))

}


# Dual Power Distortion - Speacial Case

dp_mean <- function(b,rho){

    return (beta_mean(1,b,rho))

}

dp_variance <- function(b,rho){

    return(beta_variance(1,b,rho))

}

dp_default_prob <- function(b,PD,rho){

    return(beta_default_prob(1, b, PD, rho))

}

dp_default_corr <- function(b,PD,rho){

    return(beta_default_corr(1, b, PD, rho))

}

find_rho_dp <- function(target_DDC, b, PD, tol = 1e-10){

    return(find_rho_beta(target_DDC, 1, b, PD, tol))

} 
