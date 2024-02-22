import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
source('~/dev/pabm-grdpg/functions.R')
source('~/dev/manifold-block-models/functions.R')
library(ggplot2)

doMC::registerDoMC(parallel::detectCores())

nrep <- 50

K <- 2
beta0 <- 1
beta1 <- 2
beta2 <- -1
sigma_eps <- 1

N <- 2 ^ 6
n.vec <- 2 ^ c(10, 9, 8, 7, 6)

sim.dir <- '~/dev/multilayer-rdpg/simulations/beta-regression'

parallelize.outer <- TRUE
if (parallelize.outer) {
  `%DO%` <- `%dopar%`
  parallelize.curvefit <- FALSE
} else {
  `%DO%` <- `%do%`
  parallelize.curvefit <- TRUE
}

out.df <- foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
  print(paste('n =', n))
  foreach(rep = seq(nrep), .combine = dplyr::bind_rows) %do% {
    print(paste('rep', rep, 'out of', nrep))
    rep.filename <- paste0('beta-reg-n-', n, '-rep-', rep, '.csv')
    if (rep.filename %in% dir(sim.dir)) {
      one.rep.df <- readr::read_csv(file.path(sim.dir, rep.filename),
                                    progress = FALSE,
                                    show_col_types = FALSE)
    } else {
      datamat.df <- plyr::ldply(seq(N), function(i) {
        set.seed(i + (rep - 1) * N)
        print(paste(i, '/', N))
        iter.filename <- paste0('beta-reg-n-', n, '-rep-', rep, '-iter-', i, '.csv')
        if (iter.filename %in% dir(sim.dir)) {
          one.iter.df <- readr::read_csv(file.path(sim.dir, iter.filename),
                                         progress = FALSE,
                                         show_col_types = FALSE)
        } else {
          a <- runif(1, 1, 2)
          b <- runif(1, 1, 2)
          y <- beta0 + beta1 * a + beta2 * b + rnorm(1, sigma_eps)
          t. <- rbeta(n, a, b)
          x1 <- t. ^ 2
          x2 <- 2 * t. * (1 - t.)
          X <- cbind(x1, x2)
          P <- X %*% t(X)
          A <- draw.graph(P)
          Xhat <- embedding(A, 2, 0)
          curve.est <- estimate.bezier.curve.2(Xhat, 
                                               degree = 2, 
                                               intercept = FALSE, 
                                               initialization = 'isomap', 
                                               # min.t = 0, max.t = 1,
                                               normalize = TRUE,
                                               parallel = parallelize.curvefit)
          t.hat <- curve.est$t
          param.est <- EnvStats::ebeta(t.hat)
          one.iter.df <- dplyr::tibble(y = y, 
                                       a = a, 
                                       b = b,
                                       a.hat = param.est$parameters[1], 
                                       b.hat = param.est$parameters[2])
          readr::write_csv(one.iter.df,
                           file.path(sim.dir, iter.filename))
        }
        return(one.iter.df)
      }, .parallel = parallelize.outer)
      
      # datamat.df <- foreach(i = seq(N), .combine = dplyr::bind_rows, .errorhandling = 'remove') %DO% {
      #   set.seed(i)
      #   print(paste(i, '/', N))
      #   iter.filename <- paste0('beta-reg-n-', n, '-rep-', rep, '-iter-', i, '.csv')
      #   if (iter.filename %in% dir(sim.dir)) {
      #     one.iter.df <- readr::read_csv(file.path(sim.dir, iter.filename),
      #                                    progress = FALSE,
      #                                    show_col_types = FALSE)
      #   } else {
      #   a <- runif(1, 0, 2)
      #   b <- runif(1, 0, 2)
      #   y <- beta0 + beta1 * a + beta2 * b + rnorm(1, sigma_eps)
      #   t. <- rbeta(n, a, b)
      #   x1 <- t. ^ 2
      #   x2 <- 2 * t. * (1 - t.)
      #   X <- cbind(x1, x2)
      #   P <- X %*% t(X)
      #   A <- draw.graph(P)
      #   Xhat <- embedding(A, 2, 0)
      #   curve.est <- estimate.bezier.curve.2(Xhat, 
      #                                        degree = 2, 
      #                                        intercept = FALSE, 
      #                                        initialization = 'isomap', 
      #                                        # min.t = 0, max.t = 1,
      #                                        normalize = TRUE,
      #                                        parallel = parallelize.curvefit)
      #   plot(Xhat); points(curve.est$X, col = 2)
      #   t.hat <- curve.est$t
      #   plot(t.hat, t.)
      #   param.est <- EnvStats::ebeta(t.hat)
      #   one.iter.df <- dplyr::tibble(y = y, 
      #                                a = a, 
      #                                b = b,
      #                                a.hat = param.est$parameters[1], 
      #                                b.hat = param.est$parameters[2])
      #   readr::write_csv(one.iter.df,
      #                    file.path(sim.dir, iter.filename))
      #   }
      #   return(one.iter.df)
      # }
      
      train.ind <- sample(seq(N), N / 2)
      train.df <- datamat.df[train.ind, ]
      test.df <- datamat.df[-train.ind, ]
      y.test <- test.df$y
      lm.out <- lm(y ~ a.hat + b.hat, data = train.df)
      y.hat <- predict(lm.out, test.df)
      mse <- mean((y.test - y.hat) ^ 2)
      a.mse <- mean((datamat.df$a - datamat.df$a.hat) ^ 2)
      b.mse <- mean((datamat.df$b - datamat.df$b.hat) ^ 2)
      one.rep.df <- dplyr::tibble(n = n, 
                                  mse = mse,
                                  a.mse = a.mse,
                                  b.mse = b.mse)
      readr::write_csv(one.rep.df, 
                       file.path(sim.dir, rep.filename))
    }
    return(one.rep.df)
  }
}

ggplot(out.df) + 
  geom_boxplot(aes(x = n, y = mse, group = n)) + 
  scale_x_log10(breaks = sort(n.vec)) +
  scale_y_log10() +
  labs(x = 'number of vertices', y = 'MSE')

ggplot(out.df) + 
  geom_boxplot(aes(x = n, y = mse, group = n),
               outlier.shape = NA) + 
  scale_x_log10(breaks = sort(n.vec)) +
  scale_y_log10(limits = quantile(out.df$mse, c(.1, .9))) +
  labs(x = 'number of vertices', y = 'MSE')

ggplot(out.df) + 
  geom_boxplot(aes(x = n, y = mse, group = n),
               outlier.shape = NA) + 
  scale_x_log10(breaks = sort(n.vec)) +
  scale_y_continuous(limits = quantile(out.df$mse, c(.1, .9))) +
  labs(x = 'number of vertices', y = 'MSE')

readr::write_csv(out.df, '~/dev/multilayer-rdpg/simulations/beta-regression/beta-regression-sim.csv')
