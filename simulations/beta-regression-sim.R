import::from(magrittr, `%>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
source('~/dev/pabm-grdpg/functions.R')
source('~/dev/manifold-block-models/functions.R')
library(ggplot2)

doMC::registerDoMC(parallel::detectCores())

nrep <- 50

beta0 <- 1
beta1 <- 2
beta2 <- -1
sigma_eps <- .5

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
    datamat.df <- plyr::ldply(seq(N), function(i) {
      set.seed(i + (rep - 1) * N)
      print(paste(i, '/', N))
      a <- runif(1, 1, 2)
      b <- runif(1, 1, 2)
      y <- beta0 + beta1 * a + beta2 * b + rnorm(1, 0, sigma_eps)
      t. <- rbeta(n, a, b)
      x1 <- t. ^ 2
      x2 <- 2 * t. * (1 - t.)
      X <- cbind(x1, x2)
      P <- X %*% t(X)
      A <- draw.graph(P)
      Xhat <- embedding(A, 2, 0)
      
      p.init <- matrix(c(0, 0, 
                         Xhat[which.max(apply(Xhat, 1, function(x) sum(x ^ 2))), ]),
                       ncol = 2, 
                       byrow = TRUE)
      curve.est <- estimate.bezier.curve.2(Xhat,
                                           degree = 2,
                                           intercept = FALSE,
                                           # initialization = 'isomap',
                                           init.params = p.init,
                                           # min.t = 0, max.t = 1,
                                           normalize = TRUE,
                                           parallel = parallelize.curvefit)
      t.hat <- curve.est$t
      param.est <- EnvStats::ebeta(t.hat)
      
      avg.degree <- mean(colSums(A))
      
      A.graph <- igraph::graph_from_adjacency_matrix(A)
      transitivity <- igraph::transitivity(A.graph)
      
      one.iter.df <- dplyr::tibble(n = n,
                                   y = y,
                                   a = a,
                                   b = b,
                                   a.hat = param.est$parameters[1],
                                   b.hat = param.est$parameters[2],
                                   avg.degree = avg.degree,
                                   transitivity = transitivity)
      return(one.iter.df)
    }, .parallel = parallelize.outer)
    
    train.ind <- sample(seq(N), N / 2)
    train.df <- datamat.df[train.ind, ]
    test.df <- datamat.df[-train.ind, ]
    y.test <- test.df$y
    lm.out <- lm(y ~ a.hat + b.hat, data = train.df)
    lm.d <- lm(y ~ avg.degree, data = train.df)
    lm.t <- lm(y ~ transitivity, data = train.df)
    y.hat <- predict(lm.out, test.df)
    y.hat.d <- predict(lm.d, test.df)
    y.hat.t <- predict(lm.t, test.df)
    mse <- mean((y.test - y.hat) ^ 2)
    mse.d <- mean((y.test - y.hat.d) ^ 2)
    mse.t <- mean((y.test - y.hat.t) ^ 2)
    a.mse <- mean((datamat.df$a - datamat.df$a.hat) ^ 2)
    b.mse <- mean((datamat.df$b - datamat.df$b.hat) ^ 2)
    one.rep.df <- dplyr::tibble(n = n, 
                                mse = mse,
                                a.mse = a.mse,
                                b.mse = b.mse,
                                mse.d,
                                mse.t)
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

out.df %>% 
  dplyr::group_by(n) %>% 
  dplyr::summarise(a.se = sd(a.mse) / sqrt(dplyr::n()),
                   b.se = sd(b.mse) / sqrt(dplyr::n()),
                   a.mse = mean(a.mse),
                   b.mse = mean(b.mse)) %>% 
  ggplot() + 
  # scale_y_log10() + 
  scale_x_log10() + 
  geom_line(aes(x = n, y = a.mse, colour = 'a')) + 
  geom_errorbar(aes(x = n, ymin = a.mse - a.se, ymax = a.mse + a.se, colour = 'a')) + 
  geom_line(aes(x = n, y = b.mse, colour = 'b')) + 
  geom_errorbar(aes(x = n, ymin = b.mse - b.se, ymax = b.mse + b.se, colour = 'b'))

out.df %>% 
  dplyr::group_by(n) %>% 
  dplyr::summarise(mse.beta = mean(mse),
                   se.beta = sd(mse) / sqrt(dplyr::n()),
                   se.d = sd(mse.d) / sqrt(dplyr::n()),
                   mse.d = mean(mse.d),
                   se.t = sd(mse.t) / sqrt(dplyr::n()),
                   mse.t = mean(mse.t)) %>% 
  ggplot() + 
  # scale_y_log10() + 
  scale_x_log10() + 
  geom_line(aes(x = n, y = mse.beta, colour = 'beta')) + 
  geom_errorbar(aes(x = n, 
                    ymin = mse.beta - se.beta, 
                    ymax = mse.beta + se.beta,
                    colour = 'beta')) + 
  geom_line(aes(x = n, y = mse.d, colour = 'average degree')) + 
  geom_errorbar(aes(x = n,
                    ymin = mse.d - se.d,
                    ymax = mse.d + se.d,
                    colour = 'average degree')) + 
  
  geom_line(aes(x = n, y = mse.t, colour = 'transitivity')) + 
  geom_errorbar(aes(x = n,
                    ymin = mse.t - se.t,
                    ymax = mse.t + se.t,
                    colour = 'transitivity'))

readr::write_csv(out.df, '~/dev/multilayer-rdpg/simulations/beta-regression/beta-regression-sim.csv')
