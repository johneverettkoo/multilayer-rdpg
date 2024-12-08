kl.divergence <- function(x1, x2) {
  Sigma1 <- cov(x1)
  Sigma2 <- cov(x2)
  mu1 <- colMeans(x1)
  mu2 <- colMeans(x2)
  
  d <- length(mu1)
  
  .5 * (log(det(Sigma2) / det(Sigma1)) - d + 
          sum(diag(solve(Sigma2, Sigma1))) + 
          as.numeric(t(mu2 - mu1) %*% solve(Sigma2, (mu2 - mu1))))
}


n.i <- 1024

a <- 1
b <- 1

baseline.df <- plyr::ldply(seq(100), function(i) {
  t. <- rbeta(n.i, a, b)
  x1 <- t. ^ 2
  x2 <- 2 * t. * (1 - t.)
  X <- cbind(x1, x2)
  P <- X %*% t(X)
  A <- draw.graph(P)
  
  transitivity <- A %>% 
    igraph::graph_from_adjacency_matrix() %>% 
    igraph::transitivity()
  
  avg.degree <- mean(colSums(A)) / n
  
  dplyr::tibble(transitivity = transitivity, 
                avg.degree =avg.degree)
})

summary(baseline.df)

a.mod <- seq(.5, 2.5, .1)
b.mod <- seq(.5, 2.5, .1)
search.df <- expand.grid(a = a.mod, b = b.mod)

out.df <- plyr::ldply(seq(nrow(search.df)), function(i) {
  a.mod<- search.df$a[i]
  b.mod <- search.df$b[i]
  
  mod.df <- plyr::ldply(seq(100), function(i) {
    t. <- rbeta(n.i, a.mod, b.mod)
    x1 <- t. ^ 2
    x2 <- 2 * t. * (1 - t.)
    X <- cbind(x1, x2)
    P <- X %*% t(X)
    A <- draw.graph(P)
    
    transitivity <- A %>% 
      igraph::graph_from_adjacency_matrix() %>% 
      igraph::transitivity()
    
    avg.degree <- mean(colSums(A)) / n
    
    dplyr::tibble(transitivity = transitivity, 
                  avg.degree =avg.degree)
  })
  
  train.df <- dplyr::bind_rows(baseline.df %>% dplyr::mutate(z = 1),
                               mod.df %>% dplyr::mutate(z = 2))
  model <- MASS::lda(z ~ avg.degree + transitivity, data = train.df)
  zhat <- predict(model, train.df)$posterior %>% 
    apply(1, which.max)
  acc <- mean(zhat == train.df$z)
  
  dplyr::tibble(a = a.mod,
                b = b.mod, transitivity = mean(mod.df$transitivity),
                avg.degree = mean(mod.df$avg.degree),
                kl.divergence = kl.divergence(baseline.df, mod.df),
                error = 1 - acc)
}, .parallel = TRUE)

out.df %>% 
  dplyr::mutate(dt = abs(transitivity - mean(baseline.df$transitivity)), 
                dd = abs(avg.degree - mean(baseline.df$avg.degree))) %>% 
  View()

a.mod <- 1.2
b.mod <- 1.2

mod.df <- plyr::ldply(seq(100), function(i) {
  t. <- rbeta(n.i, a.mod, b.mod)
  x1 <- t. ^ 2
  x2 <- 2 * t. * (1 - t.)
  X <- cbind(x1, x2)
  P <- X %*% t(X)
  A <- draw.graph(P)
  
  transitivity <- A %>% 
    igraph::graph_from_adjacency_matrix() %>% 
    igraph::transitivity()
  
  avg.degree <- mean(colSums(A)) / n
  
  dplyr::tibble(transitivity = transitivity, 
                avg.degree =avg.degree)
})

plot(density(baseline.df$transitivity))
lines(density(mod.df$transitivity), col = 2)

plot(density(baseline.df$avg.degree))
lines(density(mod.df$avg.degree), col = 2)
