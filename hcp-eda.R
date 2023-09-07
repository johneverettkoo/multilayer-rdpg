import::from(magrittr, `%>%`, `%<>%`)
import::from(foreach, `%do%`, `%dopar%`)
source('~/dev/pabm-grdpg/functions.R')
source('~/dev/manifold-block-models/functions.R')
library(mclust)
library(ggplot2)

doMC::registerDoMC(parallel::detectCores())

data.dir <- '~/Downloads/repeated_10_scale_125/'
graph.list <- lapply(dir(data.dir), function(x) {
  graph <- igraph::read.graph(file.path(data.dir, x), format = 'graphml')
  A <- graph %>% 
    igraph::as_adjacency_matrix() %>% 
    as.matrix()
  z <- graph %>% 
    igraph::get.vertex.attribute() %>% 
    {.$dn_hemisphere} %>% 
    as.factor()
  return(list(A = A, z = z, id = substr(x, 1, 6)))
})

datamat.df <- foreach::foreach(i = seq_along(graph.list), .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
  # print(paste(i, 'out of', length(graph.list)))
  graph <- graph.list[[i]]
  A <- graph$A
  z <- graph$z
  Xhat <- embedding(A, 2, 0, scale = TRUE)
  clustering <- manifold.clustering(Xhat, 
                                    degree = 1,
                                    initialization = as.numeric(z),
                                    intercept = FALSE)
  params.l <- EnvStats::ebeta(clustering$t[[1]])$parameters
  params.r <- EnvStats::ebeta(clustering$t[[2]])$parameters
  angle.l <- atan2(clustering$p[[1]][2, 2],
                   clustering$p[[1]][2, 1])
  angle.r <- atan2(clustering$p[[2]][2, 2],
                   clustering$p[[2]][2, 1])
  angle.diff <- abs(angle.r - angle.l)
  if (angle.diff > pi) {
    angle.diff <- 2 * pi - angle.diff
  }
  dplyr::tibble(id = graph$id,
                a.l = params.l[1],
                b.l = params.l[2],
                a.r = params.r[1],
                b.r = params.r[2],
                angle.diff = angle.diff)
}

load('~/Downloads/HCP/HCP-Y/DTI_fMRI/Data.RData')

datamat.df %<>% dplyr::left_join(
  dat.demo %>% 
    dplyr::transmute(id = as.character(Subject),
                     gender = Gender,
                     age = Age),
  by = 'id'
)

readr::write_csv(datamat.df, '~/Documents/postdoc-first-project/hcp-eda.csv')

datamat.df %>% 
  na.omit() %>% 
  GGally::ggpairs(columns = 2:6, aes(colour = gender))

datamat.df %>% 
  na.omit() %>% 
  dplyr::filter(age != '36+') %>% 
  GGally::ggpairs(columns = 2:6, aes(colour = age))
