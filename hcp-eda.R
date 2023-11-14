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

extract.t.norm <- function(roi, A, z, p, t) {
  roi.l <- paste(roi, 'L', sep = '_')
  roi.r <- paste(roi, 'R', sep = '_')
  
  idx.l <- which(colnames(A)[z == 1] == roi.l)
  idx.r <- which(colnames(A)[z == 2] == roi.r)
  
  t.l <- t[[1]][idx.l]
  t.r <- t[[2]][idx.r]
  
  norm.l <- (p[[1]] * t.l) %>% 
    magrittr::raise_to_power(2) %>% 
    sum() %>% 
    sqrt()
  norm.r <- (p[[2]] * t.r) %>% 
    magrittr::raise_to_power(2) %>% 
    sum() %>% 
    sqrt()
  
  return(list(t = list(l = t.l,
                       r = t.r),
              norm = list(l = norm.l,
                          r = norm.r)))
}

# data.dir <- '~/Downloads/repeated_10_scale_125/'
# graph.list <- lapply(dir(data.dir), function(x) {
#   graph <- igraph::read.graph(file.path(data.dir, x), format = 'graphml')
#   A <- graph %>% 
#     igraph::as_adjacency_matrix() %>% 
#     as.matrix()
#   z <- graph %>% 
#     igraph::get.vertex.attribute() %>% 
#     {.$dn_hemisphere} %>% 
#     as.factor()
#   return(list(A = A, z = z, id = substr(x, 1, 6)))
# })

load('~/Downloads/HCP/HCP-Y/DTI_fMRI/Data.RData')

# graph.list <- lapply(graph.list, function(graph) {
#   if (graph$id %in% as.character(dat.demo$Subject)) {
#     return(graph)
#   } else {
#     return(NULL)
#   }
# })
# graph.list <- graph.list[!sapply(graph.list, is.null)]
# 
# datamat.df <- foreach(i = seq_along(graph.list), .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
#   # print(paste(i, 'out of', length(graph.list)))
#   graph <- graph.list[[i]]
#   A <- graph$A
#   z <- graph$z
#   
#   Xhat <- embedding(A, 2, 0, scale = TRUE)
#   clustering <- manifold.clustering(Xhat, 
#                                     degree = 1,
#                                     initialization = as.numeric(z),
#                                     normalize = TRUE,
#                                     intercept = FALSE,
#                                     parallel = FALSE)
#   
#   params.l <- EnvStats::ebeta(clustering$t[[1]])$parameters
#   params.r <- EnvStats::ebeta(clustering$t[[2]])$parameters
#   
#   angle.l <- atan2(clustering$p[[1]][2, 2],
#                    clustering$p[[1]][2, 1])
#   angle.r <- atan2(clustering$p[[2]][2, 2],
#                    clustering$p[[2]][2, 1])
#   angle.diff <- abs(angle.r - angle.l)
#   if (angle.diff > pi) {
#     angle.diff <- 2 * pi - angle.diff
#   }
#   
#   radius.l <- clustering$p[[1]] %>% 
#     magrittr::raise_to_power(2) %>% 
#     sum() %>% 
#     sqrt()
#   radius.r <- clustering$p[[2]] %>% 
#     magrittr::raise_to_power(2) %>% 
#     sum() %>% 
#     sqrt()
#   
#   dplyr::tibble(id = graph$id,
#                 a.l = params.l[1],
#                 b.l = params.l[2],
#                 a.r = params.r[1],
#                 b.r = params.r[2],
#                 angle.diff = angle.diff,
#                 radius.l = radius.l,
#                 radius.r = radius.r)
# }

rois <- colnames(dti.countV2[, , 1])
rois.without.hemisphere <- rois %>% 
  stringr::str_sub(end = -3) %>% 
  unique()

ids <- dimnames(dti.countV2)[[3]]


datamat.df <- foreach(id = ids, .combine = dplyr::bind_rows, .errorhandling = 'remove') %dopar% {
  A <- dti.countV2[, , id] %>% 
    unique() %>% 
    t() %>% 
    unique()
  A <- A / mean(A)
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
  
  plot.estimated.curves(Xhat, clustering) +
    scale_colour_brewer(palette = 'Set1') +
    ggrepel::geom_label_repel(aes(x = Xhat[, 1],
                                  y = Xhat[, 2],
                                  label = colnames(A))) +
    ggtitle(paste0('id=', id))
  
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
  
  individual.stats.df <- dplyr::tibble(
    id = id,
    a.l = params.l[1],
    b.l = params.l[2],
    a.r = params.r[1],
    b.r = params.r[2],
    angle.diff = angle.diff,
    radius.l = radius.l,
    radius.r = radius.r
  )
  
  for (roi in rois.without.hemisphere) {
    roi.stats <- extract.t.norm(roi = roi,
                                A = A,
                                z = clustering$z,
                                p = clustering$p,
                                t = clustering$t)
    individual.stats.df[[paste0(roi, '.t.l')]] <- roi.stats$t$l
    individual.stats.df[[paste0(roi, '.t.r')]] <- roi.stats$t$r
    individual.stats.df[[paste0(roi, '.norm.l')]] <- roi.stats$norm$l
    individual.stats.df[[paste0(roi, '.norm.r')]] <- roi.stats$norm$r
    
  }
  
  return(individual.stats.df)
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
                     agreeableness = NEOFAC_A,
                     openness = NEOFAC_O,
                     conscientiousness = NEOFAC_C,
                     neuroticism = NEOFAC_N,
                     extraversion = NEOFAC_E,
                     gender = Gender,
                     age = Age),
  by = 'id'
)

readr::write_csv(datamat.df, '~/Documents/postdoc-first-project/hcp-eda.csv')

# i = 12
# 
# datamat.df %>% 
#   na.omit() %>% 
#   GGally::ggpairs(columns = c(2:8, i), aes(colour = gender))
# 
# datamat.df %>% 
#   na.omit() %>% 
#   dplyr::filter(age != '36+') %>% 
#   GGally::ggpairs(columns = 2:10, aes(colour = age))


plot(strength ~ radius.r, data = datamat.df)

responses <- c('picseq',
               'picvocab',
               'odor',
               'cardsort',
               'flanker',
               'read',
               'procspeed',
               'listsort',
               'endurance',
               'dexterity',
               'strength',
               'taste',
               'agreeableness',
               'openness',
               'conscientiousness',
               'neuroticism',
               'extraversion')
inputs <- rois.without.hemisphere %>% 
  expand.grid(c('t.l', 't.r', 'norm.l', 'norm.r')) %>% 
  apply(1, paste, collapse = '.') %>% 
  c('a.l', 'b.l', 'a.r', 'b.r', 'radius.l', 'radius.r', 'angle.diff')

cor.df <- expand.grid(x = inputs, y = responses)

cor.df$cor <- apply(cor.df, 1, function(row) {
  cor(datamat.df[[row[1]]], 
      datamat.df[[row[2]]], 
      use = 'pairwise.complete.obs')
})
cor.df %>% 
  dplyr::arrange(y, -abs(cor)) %>% 
  View()

gender.tval.df <- plyr::ldply(inputs, function(x) {
  t.val <- datamat.df %>% 
    dplyr::select(all_of(c(x, 'gender'))) %>% 
    lm() %>% 
    summary() %>% 
    magrittr::extract2('coefficients') %>% 
    {.['genderM', 't value']}
  x.scaled <- datamat.df[[x]] %>% 
    scale()
  f.median <- x.scaled[datamat.df$gender == 'F'] %>% 
    median()
  m.median <- x.scaled[datamat.df$gender == 'M'] %>% 
    median()
  
  dplyr::tibble(x = x, t.val = t.val, diff = m.median - f.median) %>% 
    return()
})

set.seed(1)
train.ind <- sample(c(TRUE, FALSE), nrow(datamat.df), replace = TRUE)
train.df <- datamat.df[train.ind, ]
test.df <- datamat.df[!train.ind, ]

y.var <- 'picvocab'

mod <- lm(formula(paste(y.var, '~ 1')), data = train.df %>% 
     dplyr::mutate_at(inputs, ~(scale(.) %>% as.vector()))) %>% 
  step(direction = 'both', 
       # k = log(nrow(train.df)),
       scope = formula(paste('~', '(', paste(inputs, collapse = '+'), ') * gender')))
summary(mod)
# plot(mod)
# yhat <- mod$fitted.values
yhat <- predict(mod, newdata = test.df)
plot(yhat, test.df[[y.var]], asp = 1, col = factor(test.df$gender)); abline(0, 1, col = 4)
cor(yhat, test.df[[y.var]], use = 'complete.obs')

# inputs <- c('radius.l', 'radius.r',
#             'a.l', 'a.r',
#             'b.l', 'b.r',
#             'angle.diff')

x.train <- train.df %>% 
  dplyr::select(inputs) %>% 
  as.matrix() %>% 
  scale() %>% 
  cbind(gender = as.numeric(as.factor(train.df$gender)) - 1)
x.test <- test.df %>% 
  dplyr::select(inputs) %>% 
  as.matrix() %>% 
  scale() %>% 
  cbind(gender = as.numeric(as.factor(test.df$gender)) - 1)
y.train <- train.df[[y.var]]
y.test <- test.df[[y.var]]
dtrain <- xgboost::xgb.DMatrix(x.train, label = y.train)
dtest <- xgboost::xgb.DMatrix(x.test, label = y.test)
xgb.mod <- xgboost::xgb.train(data = dtrain,
                              nrounds = 1000,
                              params = list(nthread = 10,
                                            # objective = 'binary:logistic',
                                            # eval_metric = 'auc',
                                            # booster = 'gblinear',
                                            max_depth = 2,
                                            subsample = .5,
                                            alpha = 1,
                                            eta = .1),
                              watchlist = list(train = dtrain,
                                               eval = dtest))
xgb.mod$evaluation_log %>% 
  dplyr::mutate(train_r2 = 1 - train_rmse ^ 2 / var(train.df[[y.var]], na.rm = TRUE),
                test_r2 = 1 - eval_rmse ^ 2 / var(test.df[[y.var]], na.rm = TRUE)) %>% 
  dplyr::filter(train_r2 > 0) %>% 
  ggplot() + 
  # scale_y_log10() + 
  geom_line(aes(x = iter, y = train_r2)) + 
  geom_line(aes(x = iter, y = test_r2), colour = 'red')

xgb.mod$evaluation_log %>% 
  ggplot() + 
  geom_line(aes(x = iter, y = train_auc)) + 
  geom_line(aes(x = iter, y = eval_auc), colour = 'red') + 
  ylim(.5, 1)

xgboost::xgb.importance(model = xgb.mod) %>% 
  xgboost::xgb.plot.importance(top_n = 10)

xgb.mod$evaluation_log %>% 
  dplyr::filter(train_rmse < sd(datamat.df[[y.var]])) %>% 
  ggplot() + 
  # scale_y_log10() + 
  geom_line(aes(x = iter, y = train_rmse)) + 
  geom_line(aes(x = iter, y = eval_rmse), colour = 'red') + 
  geom_hline(aes(yintercept = sd(datamat.df[[y.var]])), colour = 'blue')
