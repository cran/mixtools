regmixEM.mixed = function (y, x, w = NULL, sigma = NULL, arb.sigma = TRUE, alpha = NULL, 
    lambda = NULL, mu = NULL, R = NULL, arb.R = TRUE, k = 2, 
    mixed = FALSE, addintercept.fixed = FALSE, addintercept.random = TRUE, 
    epsilon = 1e-08, maxit = 10000, verb = FALSE) 
{
    `%L*%` = function(x, y) lapply(1:length(x), function(z) {
        t(x[[z]]) %*% y[[z]]
    })
    N <- length(y)
    n <- sapply(y, length)
    if (addintercept.random) {
        x.1 = lapply(1:N, function(i) cbind(1, x[[i]]))
        x = x.1
    }
    p <- ncol(x[[1]])
    if (mixed == TRUE) {
        if (addintercept.fixed) {
            w.1 = lapply(1:N, function(i) cbind(1, w[[i]]))
            w = w.1
        }
        tww = w %L*% w
        tww.inv = 0
        for (i in 1:N) {
            tww.inv = tww.inv + tww[[i]]
        }
        tww.inv = solve(tww.inv)
        twx = w %L*% x
        twy = w %L*% y
        q <- ncol(w[[1]])
    }
    if (is.null(w)) {
        w = as.list(rep(0, N))
    }
    tmp <- regmix.mixed.init(y = y, x = x, w = w, sigma = sigma, 
        arb.sigma = arb.sigma, alpha = alpha, lambda = lambda, 
        mu = mu, R = R, arb.R = arb.R, k = k, mixed = mixed, 
        addintercept.fixed = addintercept.fixed, addintercept.random = addintercept.random)
    alpha <- tmp$alpha
    lambda <- tmp$lambda
    mu <- tmp$mu
    R <- tmp$R
    sigma <- tmp$sigma
    s.2 <- sigma^2
    txx = x %L*% x
    diff <- 1
    iter <- 0
    L = matrix(nrow = N, ncol = k)
    if (arb.R) {
        L = t(sapply(1:N, function(i) {
            sapply(1:k, function(j) {
                dmvnorm(as.vector(y[[i]]), as.vector(x[[i]] %*% 
                  mu[, j] + as.vector(w[[i]] %*% alpha)), x[[i]] %*% 
                  R[[j]] %*% t(x[[i]]) + s.2[j * arb.sigma + 
                  (1 - arb.sigma)] * diag(1, n[i]))
            })
        }))
    }
    else {
        L = t(sapply(1:N, function(i) {
            sapply(1:k, function(j) {
                dmvnorm(as.vector(y[[i]]), as.vector(x[[i]] %*% 
                  mu[, j] + as.vector(w[[i]] %*% alpha)), x[[i]] %*% 
                  R %*% t(x[[i]]) + s.2[j * arb.sigma + (1 - 
                  arb.sigma)] * diag(1, n[i]))
            })
        }))
    }
    L.n = t(apply(t(L), 2, "*", matrix(lambda, nrow = 1)))
    obsloglik = sum(log(apply(L.n, 1, sum)))
    ll <- obsloglik
    while (diff > epsilon && iter < maxit) {
        if (arb.R) {
            V.beta = list()
            for (i in 1:N) {
                V.beta[[i]] = list()
                V.beta[[i]] = lapply(1:k, function(j) {
                  solve(txx[[i]]/s.2[j * arb.sigma + (1 - arb.sigma)] + 
                    solve(R[[j]]))
                })
            }
            beta = list()
            for (i in 1:N) {
                beta[[i]] = matrix(nrow = p, ncol = k)
                for (j in 1:k) {
                  beta[[i]][, j] = V.beta[[i]][[j]] %*% (t(x[[i]])/s.2[j * 
                    arb.sigma + (1 - arb.sigma)]) %*% (y[[i]] - 
                    as.vector(w[[i]] %*% alpha) - x[[i]] %*% 
                    mu[, j]) + mu[, j]
                }
            }
            z = matrix(nrow = N, ncol = k)
            for (i in 1:N) {
                for (j in 1:k) {
                  z.denom = c()
                  for (m in 1:k) {
                    z.denom = c(z.denom, lambda[m]/lambda[j] * 
                      (det(x[[i]] %*% R[[j]] %*% t(x[[i]]) + 
                        s.2[j * arb.sigma + (1 - arb.sigma)] * 
                          diag(1, n[i]))/det(x[[i]] %*% R[[m]] %*% 
                        t(x[[i]]) + s.2[m * arb.sigma + (1 - 
                        arb.sigma)] * diag(1, n[i])))^(0.5) * 
                      exp(-0.5 * (t(y[[i]] - as.vector(w[[i]] %*% 
                        alpha) - x[[i]] %*% mu[, m]) %*% solve(x[[i]] %*% 
                        R[[m]] %*% t(x[[i]]) + s.2[m * arb.sigma + 
                        (1 - arb.sigma)] * diag(1, n[i])) %*% 
                        (y[[i]] - as.vector(w[[i]] %*% alpha) - 
                          x[[i]] %*% mu[, m]) - t(y[[i]] - as.vector(w[[i]] %*% 
                        alpha) - x[[i]] %*% mu[, j]) %*% solve(x[[i]] %*% 
                        R[[j]] %*% t(x[[i]]) + s.2[j * arb.sigma + 
                        (1 - arb.sigma)] * diag(1, n[i])) %*% 
                        (y[[i]] - as.vector(w[[i]] %*% alpha) - 
                          x[[i]] %*% mu[, j]))))
                  }
                  z[i, j] = 1/sum(z.denom)
                }
            }
            sing <- sum(is.nan(z))
            sum.z = apply(z, 2, sum)
            lambda.new <- sum.z/N
            if (sum(lambda.new < 1e-08)>0) {
                sing <- 1
            }
            else {
                mu.new <- matrix(nrow = p, ncol = k)
                for (j in 1:k) {
                  mu.2 <- matrix(sapply(1:N, function(i) {
                    beta[[i]][, j]
                  }), ncol = N)
                  mu.new[, j] <- apply(t(apply(t(mu.2), 2, "*", 
                    matrix(z[, j], nrow = 1))), 1, sum)
                }
                mu.new = t(apply(t(mu.new), 2, "/", matrix(sum.z, 
                  nrow = 1)))
                if (mixed == TRUE) {
                  a.vec <- c()
                  for (i in 1:N) {
                    for (j in 1:k) {
                      a.vec = cbind(a.vec, z[i, j] * (twy[[i]] - 
                        twx[[i]] %*% beta[[i]][, j]))
                    }
                  }
                  alpha.new <- tww.inv %*% apply(a.vec, 1, sum)
                  alpha = alpha.new
                }
                s.tr <- matrix(nrow = N, ncol = k)
                z.n <- matrix(nrow = N, ncol = k)
                for (i in 1:N) {
                  for (j in 1:k) {
                    s.tr[i, j] = sum(diag(z[i, j] * ((y[[i]] - 
                      as.vector(w[[i]] %*% alpha) - x[[i]] %*% 
                      beta[[i]][, j]) %*% t(y[[i]] - as.vector(w[[i]] %*% 
                      alpha) - x[[i]] %*% beta[[i]][, j]) + x[[i]] %*% 
                      V.beta[[i]][[j]] %*% t(x[[i]]))))
                    z.n[i, j] = z[i, j] * n[i]
                  }
                }
                if (arb.sigma) {
                  s.2.new <- apply(s.tr, 2, sum)/apply(z.n, 2, 
                    sum)
                }
                else s.2.new <- sum(s.tr)/sum(n)
                R.new = list()
                for (j in 1:k) {
                  r.2 <- 0
                  for (i in 1:N) {
                    r <- z[i, j] * (V.beta[[i]][[j]] + t(t(beta[[i]][, 
                      j] - mu.new[, j])) %*% (beta[[i]][, j] - 
                      mu.new[, j]))
                    r.2 <- r.2 + r
                  }
                  R.new[[j]] = r.2/sum(z[, j])
                }
                lambda = lambda.new
                mu = mu.new
                s.2 = s.2.new
                R = R.new
                L = matrix(nrow = N, ncol = k)
                L = t(sapply(1:N, function(i) {
                  sapply(1:k, function(j) {
                    dmvnorm(as.vector(y[[i]]), as.vector(x[[i]] %*% 
                      mu[, j] + as.vector(w[[i]] %*% alpha)), 
                      x[[i]] %*% R[[j]] %*% t(x[[i]]) + s.2[j * 
                        arb.sigma + (1 - arb.sigma)] * diag(1, 
                        n[i]))
                  })
                }))
                L.n = t(apply(t(L), 2, "*", matrix(lambda, nrow = 1)))
                newobsloglik = sum(log(apply(L.n, 1, sum)))
            }
        }
        else {
            R.inv = solve(R)
            V.beta = list()
            for (i in 1:N) {
                V.beta[[i]] = list()
                V.beta[[i]] = lapply(1:k, function(j) {
                  solve(txx[[i]]/s.2[j * arb.sigma + (1 - arb.sigma)] + 
                    R.inv)
                })
            }
            beta = list()
            for (i in 1:N) {
                beta[[i]] = matrix(nrow = p, ncol = k)
                for (j in 1:k) {
                  beta[[i]][, j] = V.beta[[i]][[j]] %*% (t(x[[i]])/s.2[j * 
                    arb.sigma + (1 - arb.sigma)]) %*% (y[[i]] - 
                    as.vector(w[[i]] %*% alpha) - x[[i]] %*% 
                    mu[, j]) + mu[, j]
                }
            }
            z = matrix(nrow = N, ncol = k)
            for (i in 1:N) {
                for (j in 1:k) {
                  z.denom = c()
                  for (m in 1:k) {
                    z.denom = c(z.denom, lambda[m]/lambda[j] * 
                      (det(x[[i]] %*% R %*% t(x[[i]]) + s.2[j * 
                        arb.sigma + (1 - arb.sigma)] * diag(1, 
                        n[i]))/det(x[[i]] %*% R %*% t(x[[i]]) + 
                        s.2[m * arb.sigma + (1 - arb.sigma)] * 
                          diag(1, n[i])))^(0.5) * exp(-0.5 * 
                      (t(y[[i]] - as.vector(w[[i]] %*% alpha) - 
                        x[[i]] %*% mu[, m]) %*% solve(x[[i]] %*% 
                        R %*% t(x[[i]]) + s.2[m * arb.sigma + 
                        (1 - arb.sigma)] * diag(1, n[i])) %*% 
                        (y[[i]] - as.vector(w[[i]] %*% alpha) - 
                          x[[i]] %*% mu[, m]) - t(y[[i]] - as.vector(w[[i]] %*% 
                        alpha) - x[[i]] %*% mu[, j]) %*% solve(x[[i]] %*% 
                        R %*% t(x[[i]]) + s.2[j * arb.sigma + 
                        (1 - arb.sigma)] * diag(1, n[i])) %*% 
                        (y[[i]] - as.vector(w[[i]] %*% alpha) - 
                          x[[i]] %*% mu[, j]))))
                  }
                  z[i, j] = 1/sum(z.denom)
                }
            }
            sing <- sum(is.nan(z))
            sum.z = apply(z, 2, sum)
            lambda.new <- sum.z/N
            if (sum(lambda.new < 1e-08)>0) {
                sing <- 1
            }
            else {
                mu.new <- matrix(nrow = p, ncol = k)
                for (j in 1:k) {
                  mu.2 <- matrix(sapply(1:N, function(i) {
                    beta[[i]][, j]
                  }), ncol = N)
                  mu.new[, j] <- apply(t(apply(t(mu.2), 2, "*", 
                    matrix(z[, j], nrow = 1))), 1, sum)
                }
                mu.new = t(apply(t(mu.new), 2, "/", matrix(sum.z, 
                  nrow = 1)))
                if (mixed == TRUE) {
                  a.vec <- c()
                  for (i in 1:N) {
                    for (j in 1:k) {
                      a.vec = cbind(a.vec, z[i, j] * (twy[[i]] - 
                        twx[[i]] %*% beta[[i]][, j]))
                    }
                  }
                  alpha.new <- tww.inv %*% apply(a.vec, 1, sum)
                  alpha = alpha.new
                }
                s.tr <- matrix(nrow = N, ncol = k)
                z.n <- matrix(nrow = N, ncol = k)
                for (i in 1:N) {
                  for (j in 1:k) {
                    s.tr[i, j] = sum(diag(z[i, j] * ((y[[i]] - 
                      as.vector(w[[i]] %*% alpha) - x[[i]] %*% 
                      beta[[i]][, j]) %*% t(y[[i]] - as.vector(w[[i]] %*% 
                      alpha) - x[[i]] %*% beta[[i]][, j]) + x[[i]] %*% 
                      V.beta[[i]][[j]] %*% t(x[[i]]))))
                    z.n[i, j] = z[i, j] * n[i]
                  }
                }
                if (arb.sigma) {
                  s.2.new <- apply(s.tr, 2, sum)/apply(z.n, 2, 
                    sum)
                }
                else s.2.new <- sum(s.tr)/sum(n)
                r.3 <- 0
                for (j in 1:k) {
                  r.2 <- 0
                  for (i in 1:N) {
                    r <- z[i, j] * (V.beta[[i]][[j]] + t(t(beta[[i]][, 
                      j] - mu.new[, j])) %*% (beta[[i]][, j] - 
                      mu.new[, j]))
                    r.2 <- r.2 + r
                  }
                  r.3 <- r.3 + r.2
                }
                R.new = r.3/N
                lambda = lambda.new
                mu = mu.new
                s.2 = s.2.new
                R = R.new
                L = matrix(nrow = N, ncol = k)
                L = t(sapply(1:N, function(i) {
                  sapply(1:k, function(j) {
                    dmvnorm(as.vector(y[[i]]), as.vector(x[[i]] %*% 
                      mu[, j] + as.vector(w[[i]] %*% alpha)), 
                      x[[i]] %*% R %*% t(x[[i]]) + s.2[j * arb.sigma + 
                        (1 - arb.sigma)] * diag(1, n[i]))
                  })
                }))
                L.n = t(apply(t(L), 2, "*", matrix(lambda, nrow = 1)))
                newobsloglik = sum(log(apply(L.n, 1, sum)))
            }
        }
        if (sing > 0 || abs(newobsloglik) == Inf || is.nan(newobsloglik) || 
            sum(z) != N || newobsloglik < obsloglik) {
            cat("Need new starting values due to singularity...", 
                "\n")
            tmp <- regmix.mixed.init(y = y, x = x, w = w, arb.sigma = arb.sigma, 
                arb.R = arb.R, k = k, mixed = mixed, addintercept.fixed = addintercept.fixed, 
                addintercept.random = addintercept.random)
            alpha <- tmp$alpha
            lambda <- tmp$lambda
            mu <- tmp$mu
            R <- tmp$R
            sigma <- tmp$sigma
            s.2 <- sigma^2
            diff <- 1
            iter <- 0
            L = matrix(nrow = N, ncol = k)
            if (arb.R) {
                L = t(sapply(1:N, function(i) {
                  sapply(1:k, function(j) {
                    dmvnorm(as.vector(y[[i]]), as.vector(x[[i]] %*% 
                      mu[, j] + as.vector(w[[i]] %*% alpha)), 
                      x[[i]] %*% R[[j]] %*% t(x[[i]]) + s.2[j * 
                        arb.sigma + (1 - arb.sigma)] * diag(1, 
                        n[i]))
                  })
                }))
            }
            else {
                L = t(sapply(1:N, function(i) {
                  sapply(1:k, function(j) {
                    dmvnorm(as.vector(y[[i]]), as.vector(x[[i]] %*% 
                      mu[, j] + as.vector(w[[i]] %*% alpha)), 
                      x[[i]] %*% R %*% t(x[[i]]) + s.2[j * arb.sigma + 
                        (1 - arb.sigma)] * diag(1, n[i]))
                  })
                }))
            }
            L.n = t(apply(t(L), 2, "*", matrix(lambda, nrow = 1)))
            obsloglik = sum(log(apply(L.n, 1, sum)))
        ll <- obsloglik
        }
        else {
            diff <- newobsloglik - obsloglik
            obsloglik <- newobsloglik
        ll <- c(ll, obsloglik)
            iter <- iter + 1
            if (verb) {
                cat("iteration=", iter, "diff=", diff, "log-likelihood", 
                  obsloglik, "\n")
            }
        }
    }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, y=y, lambda = lambda, mu = mu, R = R, sigma = sqrt(s.2), alpha = alpha, 
        loglik = obsloglik, posterior.z = z, posterior.beta = beta, all.loglik = ll, ft="regmixEM.mixed")
    class(a) = "mixEM"
    a
}
