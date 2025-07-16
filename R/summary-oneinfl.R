#' Summarize One-Inflated Regression Models
#'
#' Provides a summary of the fitted model, including estimated coefficients, standard errors, 
#' significance levels, and other relevant statistics.
#'
#' @param object A model object of class `"oneinflmodel"` (for one-inflated models).
#' @param ... Additional arguments (currently unused).
#'
#' @return
#' Prints a summary table of coefficients, standard errors, z-values, p-values, 
#' significance codes, one-inflation probabilities, and log-likelihood.
#'
#' @details
#' This function generates a detailed summary of the fitted model, including:
#' - Estimated coefficients for the rate component (\code{beta}).
#' - Estimated coefficients for the one-inflation component (\code{gamma}).
#' - Standard errors.
#' - z-statistics, associated p-values, and corresponding significance codes.
#' - Average, and average absolute, one-inflation.
#' - Log-likelihood of the fitted model.
#'
#' @seealso
#' \code{\link{oneinfl}} for fitting one-inflated models.
#'
#' @examples
#' # Example usage
#' df <- data.frame(x = rnorm(100), z = rnorm(100), y = rpois(100, lambda = 5))
#' model <- oneinfl(y ~ x | z, df = df, dist = "Poisson")
#' summary(model)
#'
#' @export

summary.oneinflmodel <- function(object, ...) {
    # Extract components from the model object
    beta_vals <- object$beta
    gamma_vals <- object$gamma
    vcov_matrix_beta <- object$vc[1:length(beta_vals), 1:length(beta_vals)]
    vcov_matrix_gamma <- object$vc[(1 + length(beta_vals)):(length(beta_vals) + length(gamma_vals)), 
                                   (1 + length(beta_vals)):(length(beta_vals) + length(gamma_vals))]
    log_likelihood <- object$logl
  
    # Function to create table
    create_table <- function(coefs, vcov_matrix) {
      se <- sqrt(diag(as.matrix(vcov_matrix)))
      z_value <- coefs / se
      p_value <- 2 * (1 - pnorm(abs(z_value)))
    
      tabl <- cbind(
        Estimate = coefs,
        Std.Error = se,
        z_value = z_value,
        p.value = p_value
      )
    
      significance <- sapply(p_value, get_significance)
      result_df <- data.frame(tabl, significance)
      colnames(result_df)[5] <- ""
      return(result_df)
    }
  
    # Function to determine significance symbols
    get_significance <- function(p_value) {
      if (p_value < 0.001) {
        return("***")
      } else if (p_value < 0.01) {
        return("**")
      } else if (p_value < 0.05) {
        return("*")
      } else if (p_value < 0.1) {
        return(".")
      } else {
        return("")
      }
    }
  
    # Creating tables for beta and gamma coefficients
    beta_table <- create_table(beta_vals, vcov_matrix_beta)
    gamma_table <- create_table(gamma_vals, vcov_matrix_gamma)
  
    # Printing the results
    cat("Call:\n")
    cat(paste("formula: ", deparse(object$formula), "\n"))
    cat(paste("distribution: ", object$dist, "\n"))
  
    cat("\nCoefficients (beta):\n")
    print(beta_table, digits = 4)
  
    cat("\nCoefficients (gamma):\n")
    print(gamma_table, digits = 4)
  
    # If distribution is 'negbin', display the estimated alpha parameter
    if (object$dist == "negbin") {
      alpha_vals <- object$alpha
      vcov_matrix_alpha <- as.matrix(object$vc[nrow(object$vc), ncol(object$vc)])
      alpha_table <- create_table(alpha_vals, vcov_matrix_alpha)
      cat("\nalpha:\n")
      print(alpha_table, digits = 4)
    }
    cat(paste("\nSignif. codes:  0 `***' 0.001 `**' 0.01 `*' 0.05 `.' 0.1 ` ' 1\n"))
    cat(paste("\naverage one-inflation: ", object$avgw))
    cat(paste("\naverage absolute one-inflation: ", object$absw))
    cat(paste("\nLog-likelihood: ", log_likelihood, "\n"))
}