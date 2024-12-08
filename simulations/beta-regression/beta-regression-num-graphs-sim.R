import::from(magrittr, `%>%`, `%<>%`)
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

a.min <- .5
a.max <- 2
b.min <- .5
b.max <- 2

n <- 2 ^ 11
L.vec <- 2 ^ c(5, 6, 7)

sim.dir <- '~/dev/multilayer-rdpg/simulations/beta-regression'

parallelize.outer <- TRUE
if (parallelize.outer) {
  `%DO%` <- `%dopar%`
  parallelize.curvefit <- FALSE
} else {
  `%DO%` <- `%do%`
  parallelize.curvefit <- TRUE
}

out.df <- foreach(L = L.vec, .combine = dplyr::bind_rows) %do% {
  print(paste('L =', L))
  foreach(rep = seq(nrep), .combine = dplyr::bind_rows, .errorhandling = 'remove') %do% {
    print(paste('rep', rep, 'out of', nrep))
    datamat.df <- plyr::ldply(seq(L), function(l) {
      set.seed(l + (rep - 1) * L)
      a <- runif(1, a.min, a.max)
      b <- runif(1, a.min, a.max)
      y <- beta0 + beta1 * a + beta2 * b + rnorm(1, 0, sigma_eps)
      n.l <- rpois(1, n)
      t. <- rbeta(n.l, a, b)
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
                                           # initialization = 'x',
                                           # min.t = 0, max.t = 1,
                                           normalize = TRUE,
                                           parallel = parallelize.curvefit)
      t.hat <- curve.est$t
      param.est <- EnvStats::ebeta(t.hat)
      
      avg.degree <- mean(colSums(A)) / n.l
      
      A.graph <- igraph::graph_from_adjacency_matrix(A)
      transitivity <- igraph::transitivity(A.graph)
      
      one.iter.df <- dplyr::tibble(n = n.l,
                                   y = y,
                                   a = a,
                                   b = b,
                                   a.hat = param.est$parameters[1],
                                   b.hat = param.est$parameters[2],
                                   avg.degree = avg.degree,
                                   transitivity = transitivity)
      return(one.iter.df)
    }, .parallel = parallelize.outer)
    
    train.ind <- sample(seq(L), L / 2)
    train.df <- datamat.df[train.ind, ]
    test.df <- datamat.df[-train.ind, ]
    y.test <- test.df$y
    lm.mlsm <- lm(y ~ a.hat + b.hat, data = train.df)
    lm.stats <- lm(y ~ avg.degree + transitivity, data = train.df)
    y.hat <- predict(lm.mlsm, test.df)
    y.hat.stats <- predict(lm.stats, test.df)
    mse <- mean((y.test - y.hat) ^ 2)
    mse.stats <- mean((y.test - y.hat.stats) ^ 2)
    a.mse <- mean((datamat.df$a - datamat.df$a.hat) ^ 2)
    b.mse <- mean((datamat.df$b - datamat.df$b.hat) ^ 2)
    one.rep.df <- dplyr::tibble(n = n, 
                                L = L, 
                                mse = mse,
                                mse.stats = mse.stats,
                                a.mse = a.mse,
                                b.mse = b.mse)
    return(one.rep.df)
  }
}

trim = 0

out.df %>% 
  dplyr::group_by(L) %>% 
  dplyr::summarise(mse.beta = mean(mse, trim = trim),
                   se.beta = chemometrics::sd_trim(mse, trim = trim) / sqrt(dplyr::n()),
                   se.stats = chemometrics::sd_trim(mse.stats, trim = trim) / sqrt(dplyr::n()),
                   mse.stats = mean(mse.stats, trim = trim)) %>% 
  ggplot() + 
  scale_y_log10() +
  scale_x_log10() + 
  geom_line(aes(x = L, y = mse.beta, colour = 'mlsm')) + 
  geom_point(aes(x = L, y = mse.beta, colour = 'mlsm')) + 
  geom_errorbar(aes(x = L, 
                    ymin = mse.beta - se.beta, 
                    ymax = mse.beta + se.beta,
                    colour = 'mlsm')) + 
  geom_line(aes(x = L, y = mse.stats, colour = 'stats')) + 
  geom_point(aes(x = L, y = mse.stats, colour = 'stats')) + 
  geom_errorbar(aes(x = L,
                    ymin = mse.stats - se.stats,
                    ymax = mse.stats + se.stats,
                    colour = 'stats'))

out.df %>% 
  dplyr::group_by(L) %>% 
  dplyr::summarise(mse.mlsm.median = median(mse),
                   mse.mlsm.first.q = quantile(mse, .25),
                   mse.mlsm.third.q = quantile(mse, .75),
                   mse.stats.median = median(mse.stats),
                   mse.stats.first.q = quantile(mse.stats, .25),
                   mse.stats.third.q = quantile(mse.stats, .75)) %>% 
  ggplot() + 
  scale_y_log10() +
  scale_x_log10() + 
  geom_line(aes(x = L, y = mse.mlsm.median, colour = 'mlsm')) + 
  geom_point(aes(x = L, y = mse.mlsm.median, colour = 'mlsm')) + 
  geom_errorbar(aes(x = L, 
                    ymin = mse.mlsm.first.q, 
                    ymax = mse.mlsm.third.q,
                    colour = 'mlsm')) + 
  geom_line(aes(x = L, y = mse.stats.median, colour = 'stats')) + 
  geom_point(aes(x = L, y = mse.stats.median, colour = 'stats')) + 
  geom_errorbar(aes(x = L,
                    ymin = mse.stats.first.q,
                    ymax = mse.stats.third.q,
                    colour = 'stats'))

out.df %>% 
  dplyr::group_by(L) %>% 
  dplyr::summarise(a.se = sd(a.mse) / sqrt(dplyr::n()),
                   b.se = sd(b.mse) / sqrt(dplyr::n()),
                   a.mse = mean(a.mse),
                   b.mse = mean(b.mse)) %>% 
  ggplot() + 
  scale_y_log10() +
  scale_x_log10() + 
  geom_line(aes(x = L, y = a.mse, colour = 'a')) + 
  geom_errorbar(aes(x = L, ymin = a.mse - a.se, ymax = a.mse + a.se, colour = 'a')) + 
  geom_line(aes(x = L, y = b.mse, colour = 'b')) + 
  geom_errorbar(aes(x = L, ymin = b.mse - b.se, ymax = b.mse + b.se, colour = 'b'))

readr::write_csv(out.df, '~/dev/multilayer-rdpg/simulations/beta-regression/beta-regression-num-g
                 raphs-sim.csv')
