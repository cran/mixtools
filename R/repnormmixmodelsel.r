repnormmixmodel.sel <- function (x, k = 2, ...) 
{
    aic <- NULL
    bic <- NULL
    caic <- NULL
    icl <- NULL
        AIC <- function(emout) {
            emout$loglik - (length(emout$mu) + length(emout$stdev) + 
                length(emout$lambda) - 1)
        }
        BIC <- function(emout) {
            emout$loglik - log(nrow(x)) * (length(emout$mu) + 
                length(emout$stdev) + length(emout$lambda) - 
                1)/2
        }
        CAIC <- function(emout) {
            emout$loglik - (log(nrow(x)) + 1) * (length(emout$mu) + 
                length(emout$stdev) + length(emout$lambda) - 
                1)/2
        }
        ICL <- function(emout) {
            BIC(emout) - sum(emout$lambda * log(emout$lambda))
        }
        for (i in 1:k) {
            if (i == 1) {
                mu <- mean(x)
                loglik <- log(prod(dnorm(x, mean = mu, sd = sd(x))))
                emout <- list(mu = mu, stdev = sd(x), 
                  lambda = 1, loglik = loglik)
            }
            else emout <- repnormmixEM(x, k = i, ...)
            aic[i] <- AIC(emout)
            bic[i] <- BIC(emout)
            caic[i] <- CAIC(emout)
            icl[i] <- ICL(emout)
        }
    out = rbind(aic, bic, caic, icl)
    Winner = apply(out, 1, function(x) (1:length(x))[x == max(x)])
    rownames(out) = c("AIC", "BIC", "CAIC", "ICL")
    colnames(out) = 1:k
    cbind(out, Winner)
}
