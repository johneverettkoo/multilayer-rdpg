import::from(magrittr, `%>%`)
import::from(foreach, `%do%`, `%dopar%`)
source('~/dev/pabm-grdpg/functions.R')
source('~/dev/manifold-block-models/functions.R')
library(ggplot2)

doMC::registerDoMC(parallel::detectCores())

nrep <- 50

a.vec <- c(1, 1)
b.vec <- c(1, 1)
theta.min <- pi / 6
theta.max <- pi / 3
beta1 <- 1
beta0 <- 2
sigma_eps <- 1

N <- 2 ^ 6
n.vec <- 2 ^ c(10, 9, 8, 7, 6)

sim.dir <- '~/dev/multilayer-rdpg/simulations/angle-regression'

out.df <- foreach::foreach(n = n.vec, .combine = dplyr::bind_rows) %do% {
  print(paste('n =', n))
  foreach::foreach(rep = seq(nrep), .combine = dplyr::bind_rows) %do% {
    print(paste('rep', rep, 'out of', nrep))
    rep.filename <- paste0('angle-reg-n-', n, '-rep-', rep, '.csv')
    if (rep.filename %in% dir(sim.dir)) {
      one.rep.df <- readr::read_csv(file.path(sim.dir, rep.filename),
                                    progress = FALSE,
                                    show_col_types = FALSE)
    } else {
      datamat.df <- foreach::foreach(i = seq(N), .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
        set.seed(i + (nrep - 1) * rep)
        print(paste(i, '/', N))
        iter.filename <- paste0('angle-reg-n-', n, '-rep-', rep, '-iter-', i, '.csv')
        if (iter.filename %in% dir(sim.dir)) {
          one.iter.df <- readr::read_csv(file.path(sim.dir, iter.filename),
                                         progress = FALSE,
                                         show_col_types = FALSE)
        } else {
          z <- sample(seq(2), n, replace = TRUE)
          t. <- rbeta(n, a.vec[z], b.vec[z])
          theta <- runif(1, theta.min, theta.max)
          y <- beta0 + beta1 * theta + rnorm(1, 0, sigma_eps)
          x1 <- ifelse(z == 1, t., cos(theta) * t.)
          x2 <- ifelse(z == 1, 0, sin(theta) * t.)
          X <- cbind(x1, x2)
          P <- X %*% t(X)
          A <- draw.graph(P)
          Xhat <- embedding(A, 2, 0)
          clustering <- manifold.clustering(Xhat, 
                                            degree = 1,
                                            initialization = z,
                                            normalize = TRUE,
                                            intercept = FALSE)
          angle.1 <- atan2(clustering$p[[1]][2, 2],
                           clustering$p[[1]][2, 1])
          angle.2 <- atan2(clustering$p[[2]][2, 2],
                           clustering$p[[2]][2, 1])
          theta.hat <- abs(angle.2 - angle.1)
          if (theta.hat > pi) {
            theta.hat <- 2 * pi - theta.hat
          }
          one.iter.df <- dplyr::tibble(y = y,
                                       theta = theta,
                                       theta.hat = theta.hat)
          readr::write_csv(one.iter.df, 
                           file.path(sim.dir, iter.filename))
        }
      }
      
      train.ind <- sample(seq(N), N / 2)
      train.df <- datamat.df[train.ind, ]
      test.df <- datamat.df[-train.ind, ]
      model <- lm(y ~ theta.hat, data = train.df)
      y.test <- test.df$y
      y.hat <- predict(model, test.df)
      mse <- mean((y.test - y.hat) ^ 2)
      theta.mse <- mean((datamat.df$theta - datamat.df$theta.hat) ^ 2)
      one.rep.df <- dplyr::tibble(n = n, 
                                  mse = mse,
                                  theta.mse = theta.mse)
      readr::write_csv(one.rep.df,
                       file.path(sim.dir, rep.filename))
      return(one.rep.df)
    }
  }
}

ggplot(out.df) + 
  geom_boxplot(aes(x = n, y = mse, group = n)) + 
  scale_x_log10(breaks = sort(n.vec)) +
  scale_y_log10() +
  labs(x = 'number of vertices', y = 'MSE')

ggplot(out.df) + 
  geom_boxplot(aes(x = n, y = mse * N, group = n)) + 
  scale_x_log10(breaks = sort(n.vec)) +
  scale_y_log10() +
  labs(x = 'number of vertices', y = 'SSE')

readr::write_csv(out.df, '~/dev/multilayer-rdpg/simulations/angle-regression/angle-regression-sim.csv')
