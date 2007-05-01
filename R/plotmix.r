plot.mixEM=function (x, loglik = TRUE, density = FALSE, w = 1, alpha = 0.05, 
    marginal = FALSE, ...) 
{
    def.par <- par(no.readonly = TRUE)
    mix.object <- x
    if (!inherits(mix.object, "mixEM")) 
        stop("Use only with \"mixEM\" objects!")
    if (loglik) {
        plot(mix.object$all.loglik, col = 0, xlab = "Iteration", 
            ylab = "Log-Likelihood", main = "Observed Data Log-Likelihood", 
            ...)
        points(mix.object$all.loglik, type = "l", col = 1, lwd = 2)
    }
    if (mix.object$ft == "logisregmixEM" && density == TRUE) {
        if (ncol(mix.object$x) != 2) {
            stop(paste("The predictors must have 2 columns!", 
                "\n"))
        }
        if (sum((mix.object$y == 1) + (mix.object$y == 0)) != 
            length(mix.object$y)) {
            stop(paste("The response must be binary!", "\n"))
        }
        get(getOption("device"))()
        k = ncol(mix.object$beta)
        x = mix.object$x[, 2]
        plot(x, mix.object$y, col = (apply(mix.object$posterior, 
            1, which.max) + 1), xlab = "Predictor", ylab = "Response", 
            main = "Most Probable Component Membership")
        a = cbind(x, mix.object$y)
        a = a[order(a[, 1]), ]
        for (i in 1:k) {
            lines(a[, 1], plogis(mix.object$beta[1, i] + mix.object$beta[2, 
                i] * a[, 1]), col = (i + 1))
        }
        par(def.par)
    }
    if (mix.object$ft == "normalmixEM" && density == TRUE) {
        get(getOption("device"))()
        k <- ncol(mix.object$posterior)
        x <- sort(mix.object$x)
        a = hist(x, prob = TRUE, main = "Density Curves", plot = FALSE)
        hist(x, prob = TRUE, main = "Density Curves", ylim = c(0, 
            max(a$density) * w))
        if (length(mix.object$mu) == 1) {
            arbvar <- TRUE
            mix.object$sigma = mix.object$scale * mix.object$sigma
            arbmean <- FALSE
        }
        if (length(mix.object$mu) == k && length(mix.object$sigma) == 
            1) {
            arbmean <- TRUE
            arbvar <- FALSE
        }
        if (length(mix.object$sigma) == k && length(mix.object$mu) == 
            k) {
            arbmean <- TRUE
            arbvar <- TRUE
        }
        for (i in 1:k) {
            lines(x, mix.object$lambda[i] * dnorm(x, mean = mix.object$mu[i * 
                arbmean + (1 - arbmean)], sd = mix.object$sigma[i * 
                arbvar + (1 - arbvar)]), col = (i + 1), lwd = 2)
        }
        par(def.par)
    }
    if (mix.object$ft == "repnormmixEM" && density == TRUE) {
        get(getOption("device"))()
        x <- as.vector(as.matrix(x))
        k <- ncol(mix.object$posterior)
        x <- sort(mix.object$x)
        a = hist(x, prob = TRUE, main = "Density Curves", plot = FALSE)
        hist(x, prob = TRUE, main = "Density Curves", ylim = c(0, 
            max(a$density) * w))
        if (length(mix.object$mu) == 1) {
            arbvar <- TRUE
            mix.object$sigma = mix.object$scale * mix.object$sigma
            arbmean <- FALSE
        }
        if (length(mix.object$mu) == k && length(mix.object$sigma) == 
            1) {
            arbmean <- TRUE
            arbvar <- FALSE
        }
        if (length(mix.object$sigma) == k && length(mix.object$mu) == 
            k) {
            arbmean <- TRUE
            arbvar <- TRUE
        }
        for (i in 1:k) {
            lines(x, mix.object$lambda[i] * dnorm(x, mean = mix.object$mu[i * 
                arbmean + (1 - arbmean)], sd = mix.object$sigma[i * 
                arbvar + (1 - arbvar)]), col = (i + 1), lwd = 2)
        }
        par(def.par)
    }
    if (mix.object$ft == "regmixEM.mixed" && density == TRUE) {
        get(getOption("device"))()
        x.1 = mix.object$x
        n = sum(sapply(x.1, nrow))
        x.1.sum = sum(sapply(1:length(x.1), function(i) length(x.1[[i]][, 
            1])))
        if (x.1.sum == n) {
            x = lapply(1:length(x.1), function(i) matrix(x.1[[i]][, 
                -1], ncol = 1))
        }
        else {
            x = x.1
        }
        post.beta(x = x, y = mix.object$y, p.beta = mix.object$posterior.beta, 
            p.z = mix.object$posterior.z)
        par(def.par)
    }
    if (mix.object$ft == "mvnormalmixEM" && density == TRUE) {
        x = mix.object$x
        if (ncol(x) != 2) {
            stop(paste("The data must have 2 columns!", "\n"))
        }
        post = apply(mix.object$posterior, 1, which.max)
        k <- ncol(mix.object$posterior)
        if (is.list(mix.object$sigma)) {
            sigma = mix.object$sigma
        }
        else {
            sigma = lapply(1:k, function(i) mix.object$sigma)
        }
        if (is.list(mix.object$mu)) {
            mu = mix.object$mu
        }
        else {
            mu = lapply(1:k, function(i) mix.object$mu)
        }
        get(getOption("device"))()
        if (marginal == FALSE) {
            plot(x, col = 0)
            plot(x, col = (post + 1), xlab = "X.1", ylab = "X.2", 
                main = "Most Probable Component Membership", 
                ...)
            lapply(1:k, function(i) points(mu[[i]][1], mu[[i]][2], 
                pch = 19))
            for (i in 1:k) {
                for (j in 1:length(alpha)) {
                  ellipse(mu = mu[[i]], sigma = sigma[[i]], alpha = alpha[j], 
                    col = (i + 1))
                }
            }
        }
        else {
            x <- mix.object$x[, 1]
            y <- mix.object$x[, 2]
            xhist <- hist(x, plot = FALSE)
            yhist <- hist(y, plot = FALSE)
            top <- max(c(xhist$counts, yhist$counts))
            xrange <- range(x)
            yrange <- range(y)
            nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), 
                c(4, 1), c(1, 4), TRUE)
            layout.show(nf)
            par(mar = c(3, 3, 1, 1))
            plot(mix.object$x[, 1], col = (post + 1), mix.object$x[, 
                2], xlab = "X.1", ylab = "X.2", ...)
            lapply(1:k, function(i) points(mu[[i]][1], mu[[i]][2], 
                pch = 19))
            for (i in 1:k) {
                for (j in 1:length(alpha)) {
                  ellipse(mu = mu[[i]], sigma = sigma[[i]], alpha = alpha[j], 
                    col = (i + 1))
                }
            }
            par(mar = c(0, 3, 1, 1))
            barplot(xhist$counts, axes = FALSE, ylim = c(0, top), 
                space = 0)
            par(mar = c(3, 0, 1, 1))
            barplot(yhist$counts, axes = FALSE, xlim = c(0, top), 
                space = 0, horiz = TRUE)
        }
        par(def.par)
    }
    if (mix.object$ft == "regmixEM" && density == TRUE) {
        if (ncol(mix.object$x) != 2) {
            stop(paste("The predictors must have 2 columns!", 
                "\n"))
        }
        post = apply(mix.object$posterior, 1, which.max)
        k <- ncol(mix.object$posterior)
        x = mix.object$x[, 2]
        y = mix.object$y
        n = length(y)
        get(getOption("device"))()
        plot(x, y, col = 0, main = "Most Probable Component Membership", 
            xlab = "Predictor", ylab = "Response", ...)
        a = cbind(mix.object$x[, 2], mix.object$y, post)
        for (i in 1:k) {
            xy = subset(cbind(a, mix.object$posterior[, i]), 
                a[, 3] == i)[, -3]
            xy = matrix(xy, ncol=3)
            points(xy[, 1], xy[, 2], col = (i + 1))
            if (is.matrix(mix.object$beta) == FALSE) {
                abline(coef = mix.object$beta)
                beta = matrix(mix.object$beta, ncol = k, nrow = 2)
            }
            else {
                abline(coef = mix.object$beta[, i], col = (i + 
                  1))
		    beta = mix.object$beta
            }
            out = lm(y ~ x, weights = mix.object$posterior[, 
                i])
            fit = beta[1, i] + beta[2, i] * x
            out.aov = anova(out)
            MSE = out.aov$Mean[2]
            xy.f = cbind(x, y, fit)
            xy.sort = xy.f[order(xy.f[, 1]), ]
		x.new = seq(from=min(x),to=max(x),length.out=100)
		y.new = beta[1, i] + beta[2, i] * x.new
            s.h <- sqrt(MSE * (1/n + (x.new - mean(xy.sort[, 
                1]))^2/var(xy.sort[, 1])/(n - 1)))
            for (j in 1:length(alpha)) {
                W = sqrt(qf(1 - alpha[j], 2, n - 2))
                upper = y.new + W * s.h
                lower = y.new - W * s.h
                lines(x.new, upper, col = (i + 1))
                lines(x.new, lower, col = (i + 1))
            }
        }
        par(def.par)
    }
}