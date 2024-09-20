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

K <- 2
a.vec <- c(2, 1)
b.vec <- c(2, .839214517857)

N <- 2 ^ 6
n.vec <- 2 ^ c(11, 10, 9, 8)

sim.dir <- '~/dev/multilayer-rdpg/simulations/beta-classification'

out.df <- foreach::foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
  print(paste('n = ', n))
  foreach::foreach(rep = seq(nrep), .combine = dplyr::bind_rows) %do% {
    print(paste('rep', rep, 'out of', nrep))
    z.vec <- sample(seq(K), N, replace = TRUE)
    # a.vec <- runif(2, .25, 3)
    # b.vec <- runif(2, .25, 3)
    datamat.df <- foreach::foreach(i = seq(N), .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
      z <- z.vec[i]
      a <- a.vec[z]
      b <- b.vec[z]
      n.i <- rpois(1, n)
      t. <- rbeta(n.i, a, b)
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
                                           parallel = FALSE)
      t.hat <- curve.est$t
      param.est <- EnvStats::ebeta(t.hat)
      
      avg.degree <- mean(colSums(A))
      A.graph <- igraph::graph_from_adjacency_matrix(A)
      transitivity <- igraph::transitivity(A.graph)
      
      dplyr::tibble(z = z, 
                    a = a.vec[z],
                    b = b.vec[z],
                    a.hat = param.est$parameters[1], 
                    b.hat = param.est$parameters[2],
                    avg.degree = avg.degree,
                    transitivity = transitivity)
    }
    
    train.ind <- sample(seq(N), N / 2)
    train.df <- datamat.df[train.ind, ]
    test.df <- datamat.df[-train.ind, ]
    z.test <- z.vec[-train.ind]
    lda.mlsm <- MASS::lda(z ~ a.hat + b.hat, data = train.df)
    lda.stats <- MASS::lda(z ~ avg.degree + transitivity, data = train.df)
    z.hat.mlsm <- predict(lda.mlsm, test.df)$posterior %>% 
      apply(1, which.max)
    z.hat.stats <- predict(lda.stats, test.df)$posterior %>% 
      apply(1, which.max)
    acc.mlsm <- mean(z.test == z.hat.mlsm)
    acc.stats <- mean(z.test == z.hat.stats)
    dplyr::tibble(n = n, acc.mslm = acc.mlsm, acc.stats = acc.stats)
  }
}

ggplot(out.df) + 
  geom_boxplot(aes(x = n, y = 1 - acc.mlsm, group = n)) + 
  scale_x_log10(breaks = sort(n.vec)) + 
  scale_y_log10()

out.df %>% 
  dplyr::group_by(n) %>% 
  dplyr::summarise(acc.mslm.mean = mean(acc.mslm),
                   acc.stats.mean = mean(acc.stats),
                   acc.mslm.se = sd(acc.mslm) / dplyr::n(),
                   acc.stats.se = sd(acc.stats) / dplyr::n()) %>% 
  ggplot() + 
  scale_x_log10(breaks = sort(n.vec)) + 
  geom_line(aes(x = n, y = acc.mslm.mean, colour = 'mslm')) + 
  geom_line(aes(x = n, y = acc.stats.mean, colour = 'stats')) + 
  geom_errorbar(aes(x = n, 
                    ymin = acc.mslm.mean - 2 * acc.mslm.se, 
                    ymax = acc.mslm.mean + 2 * acc.mslm.se,
                    colour = 'mslm')) + 
  geom_errorbar(aes(x = n, 
                    ymin = acc.stats.mean - 2 * acc.stats.se, 
                    ymax = acc.stats.mean + 2 * acc.stats.se,
                    colour = 'stats'))

readr::write_csv(out.df, file.path(sim.dir, 'beta.classification-sim.csv'))
