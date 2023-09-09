import::from(magrittr, `%>%`, `%<>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
source('~/dev/pabm-grdpg/functions.R')
source('~/dev/manifold-block-models/functions.R')
library(mclust)
library(ggplot2)

doMC::registerDoMC(parallel::detectCores())

betadist.mom <- function(x) {
  xbar <- mean(x)
  s2 <- var(x)
  norm.stat <- xbar * (1 - xbar) / s2 - 1
  alpha <- xbar * norm.stat
  beta <- (1 - xbar) * norm.stat
  return(list(alpha = alpha, 
              beta = beta))
}

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

load('~/Downloads/HCP/HCP-Y/DTI_fMRI/Data.RData')

graph.list <- lapply(graph.list, function(graph) {
  if (graph$id %in% as.character(dat.demo$Subject)) {
    return(graph)
  } else {
    return(NULL)
  }
})
graph.list <- graph.list[!sapply(graph.list, is.null)]

datamat.df <- foreach(i = seq_along(graph.list), .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
  # print(paste(i, 'out of', length(graph.list)))
  graph <- graph.list[[i]]
  A <- graph$A
  z <- graph$z
  
  Xhat <- embedding(A, 2, 0, scale = TRUE)
  clustering <- manifold.clustering(Xhat, 
                                    degree = 1,
                                    initialization = as.numeric(z),
                                    normalize = TRUE,
                                    intercept = FALSE,
                                    parallel = FALSE)
  
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
  
  radius.l <- clustering$p[[1]] %>% 
    magrittr::raise_to_power(2) %>% 
    sum() %>% 
    sqrt()
  radius.r <- clustering$p[[2]] %>% 
    magrittr::raise_to_power(2) %>% 
    sum() %>% 
    sqrt()
  
  dplyr::tibble(id = graph$id,
                a.l = params.l[1],
                b.l = params.l[2],
                a.r = params.r[1],
                b.r = params.r[2],
                angle.diff = angle.diff,
                radius.l = radius.l,
                radius.r = radius.r)
}

ids <- dimnames(dti.countV2)[[3]]
datamat.df <- foreach(id = ids, .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
  A <- dti.countV2[, , id] %>% 
    unique() %>% 
    t() %>% 
    unique()
  z <- colnames(A) %>% 
    stringi::stri_sub(-1, -1) %>% 
    as.factor()
  
  Xhat <- embedding(A, 2, 0, scale = TRUE)
  clustering <- manifold.clustering(Xhat, 
                                    degree = 1,
                                    initialization = as.numeric(z),
                                    normalize = TRUE,
                                    intercept = FALSE,
                                    parallel = FALSE)
  
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
  
  radius.l <- clustering$p[[1]] %>% 
    magrittr::raise_to_power(2) %>% 
    sum() %>% 
    sqrt()
  radius.r <- clustering$p[[2]] %>% 
    magrittr::raise_to_power(2) %>% 
    sum() %>% 
    sqrt()
  
  dplyr::tibble(id = id,
                a.l = params.l[1],
                b.l = params.l[2],
                a.r = params.r[1],
                b.r = params.r[2],
                angle.diff = angle.diff,
                radius.l = radius.l,
                radius.r = radius.r)
}

datamat.df %<>% dplyr::left_join(
  dat.demo %>% 
    dplyr::transmute(id = as.character(Subject),
                     picseq = PicSeq_AgeAdj,
                     picvocab = PicVocab_AgeAdj,
                     odor = Odor_AgeAdj,
                     cardsort = CardSort_AgeAdj,
                     flanker = Flanker_AgeAdj,
                     read = ReadEng_AgeAdj,
                     procspeed = ProcSpeed_AgeAdj,
                     listsort = ListSort_AgeAdj,
                     endurance = Endurance_AgeAdj,
                     dexterity = Dexterity_AgeAdj,
                     strength = Strength_AgeAdj,
                     taste = Taste_AgeAdj,
                     gender = Gender,
                     age = Age),
  by = 'id'
)

readr::write_csv(datamat.df, '~/Documents/postdoc-first-project/hcp-eda.csv')

i = 12

datamat.df %>% 
  na.omit() %>% 
  GGally::ggpairs(columns = c(2:8, i), aes(colour = gender))

datamat.df %>% 
  na.omit() %>% 
  dplyr::filter(age != '36+') %>% 
  GGally::ggpairs(columns = 2:10, aes(colour = age))
