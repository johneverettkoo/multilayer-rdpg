import::from(magrittr, `%>%`, `%<>%`)
import::from(foreach, foreach, `%do%`, `%dopar%`)
source('~/dev/pabm-grdpg/functions.R')
source('~/dev/manifold-block-models/functions.R')
library(mclust)
library(ggplot2)

doMC::registerDoMC(parallel::detectCores())

load('~/Downloads/Data/HCP/HCP-A/DTI/Data.RData')

subject.df <- demo.data %>% 
  dplyr::filter(QC_Issue_Codes == '',
                Full_MR_Compl == 1) %>% 
  dplyr::transmute(id = src_subject_id,
                   sex = sex,
                   age = as.numeric(interview_age)) %>% 
  dplyr::filter(age < 1200) %>% 
  dplyr::mutate(age = age / 12)

z <- dti.info$hemisphere %>% 
  as.factor()

rois <- dti.info$aparc_aseg_region.names.short

A <- dat.NumFibers[[2]] %>% 
  as.matrix()
A <- (A + t(A)) / 2
diag(A) <- 0

qgraph::qgraph(A, layout = 'spring', groups = dti.info$hemisphere)
Xhat <- embedding(A, 2, 0, scale = TRUE)
clustering <- manifold.clustering(Xhat, 
                                  degree = 1,
                                  initialization = as.numeric(z),
                                  normalize = TRUE,
                                  intercept = FALSE,
                                  parallel = FALSE)
plot.estimated.curves(Xhat, clustering) + 
  ggrepel::geom_label_repel(aes(x = Xhat[, 1],
                                y = Xhat[, 2],
                                label = dti.info$aparc_aseg_region.names.short))
params.l <- EnvStats::ebeta(clustering$t[[1]])$parameters
params.r <- EnvStats::ebeta(clustering$t[[2]])$parameters

A.within.l <- A[z == 'L', z == 'L']
A.within.r <- A[z == 'R', z == 'R']
deg.within <- mean(c(A.within.l[upper.tri(A.within.l)],
                     A.within.r[upper.tri(A.within.r)]))
A.between <- A[z == 'L', z == 'R']
deg.between <- mean(A.between[upper.tri(A.between)])

ids <- intersect(names(dat.NumFibers), subject.df$id)
ase.stats.df <- foreach(id = ids, .combine = dplyr::bind_rows) %dopar% {
  A <- dat.NumFibers[[id]] %>% 
    as.matrix()
  A <- (A + t(A)) / 2
  diag(A) <- 0
  A <- A / mean(A[upper.tri(A)])
  
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
  
  A.within.l <- A[z == 'L', z == 'L']
  A.within.r <- A[z == 'R', z == 'R']
  deg.within <- mean(c(A.within.l[upper.tri(A.within.l)],
                       A.within.r[upper.tri(A.within.r)]))
  A.between <- A[z == 'L', z == 'R']
  deg.between <- mean(A.between[upper.tri(A.between)])
  
  A.eigen <- eigen(A, symmetric = TRUE)
  
  degrees <- colSums(A)
  n.l <- sum(z == 'L')
  n.r <- sum(z == 'R')
  B.ll <- 1 / n.l / n.l * sum(A[z == 'L', z == 'L'])
  B.rr <- 1 / n.r / n.r * sum(A[z == 'R', z == 'R'])
  B.lr <- 1 / n.l / n.r * sum(A[z == 'L', z == 'R'])
  
  individual.stats.df <- dplyr::tibble(
    id = id,
    # a.l = params.l[1],
    # b.l = params.l[2],
    # a.r = params.r[1],
    # b.r = params.r[2],
    angle.diff = angle.diff,
    # radius.l = radius.l,
    # radius.r = radius.r,
    deg.within = deg.within,
    deg.between = deg.between,
    lambda1 = A.eigen$values[1],
    lambda2 = A.eigen$values[2],
    B.ll = B.ll,
    B.rr = B.rr,
    B.lr = B.lr
  )
  
  return(individual.stats.df)
} %>% 
  # dplyr::filter(!(id %in% ids.to.remove)) %>% 
  dplyr::inner_join(
    subject.df, 
    by = 'id'
  )

GGally::ggpairs(ase.stats.df, columns = c(2:10), 
                aes(colour = sex)) + 
  theme_bw()

readr::write_csv(ase.stats.df,
                 '~/dev/multilayer-rdpg/hcp-a-eda.csv')

pairwise.angles.df <- foreach(i = ids, .combine = dplyr::bind_rows) %dopar% {
  A <- dat.NumFibers[[i]] %>% 
    as.matrix()
  A <- (A + t(A)) / 2
  diag(A) <- 0
  A <- A / mean(A[upper.tri(A)])
  
  Xhat <- embedding(A, 2, 0, scale = TRUE)
  
  y <- Xhat / sqrt(rowSums(Xhat ^ 2))
  
  angles <- y %>% 
    magrittr::set_rownames(colnames(A)) %>% 
    tcrossprod() %>% 
    acos()
  diag(angles) <- 0
  
  angles.df <- angles %>% 
    tibble::as_tibble() %>% 
    tibble::rownames_to_column('roi1') %>% 
    tidyr::pivot_longer(-roi1, names_to = 'roi2', values_to = 'angle') %>% 
    dplyr::mutate(roi1 = colnames(A)[as.numeric(roi1)]) %>% 
    dplyr::mutate(roi1 = factor(roi1, levels = colnames(A)),
                  roi2 = factor(roi2, levels = colnames(A))) %>% 
    dplyr::mutate(id = i,
                  age = subject.df[subject.df$id == i, ]$age)
  
  # ggplot(angles.df) + 
  #   geom_tile(aes(x = roi1, y = roi2, fill = angle)) + 
  #   viridis::scale_fill_viridis() + 
  #   coord_fixed() + 
  #   theme_minimal() + 
  #   labs(x = NULL, y = NULL) + 
  #   theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
  
  return(angles.df)
} %>% 
  dplyr::left_join(
    dti.info %>% 
      dplyr::transmute(roi1 = aparc_aseg_region.names,
                       rsn1 = RSN)
  ) %>% 
  dplyr::left_join(
    dti.info %>% 
      dplyr::transmute(roi2 = aparc_aseg_region.names,
                       rsn2 = RSN)
  ) %>% 
  dplyr::left_join(subject.df)

readr::write_csv(pairwise.angles.df, '~/dev/multilayer-rdpg/pairwise-angles.csv')

pairwise.angles.df %>% 
  dplyr::mutate(hemisphere1 = ifelse(grepl('(rh-)|(Right-)', roi1), 'R', 'L'),
                hemisphere2 = ifelse(grepl('(rh-)|(Right-)', roi2), 'R', 'L')) %>% 
  dplyr::mutate(roi1 = gsub('(ctx-rh-)|(Right-)|(ctx-lh-)|(Left-)', '', roi1),
                roi2 = gsub('(ctx-rh-)|(Right-)|(ctx-lh-)|(Left-)', '', roi2)) %>% 
  dplyr::filter(roi1 == roi2, hemisphere1 != hemisphere2) %>% 
  ggplot() + 
  geom_point(aes(x = age, y = angle, groups = roi1),
             size = .1) +
  facet_wrap(~ roi1) + 
  theme_bw()

pairwise.angles.df %>% 
  dplyr::mutate(rsn1 == rsn2) %>% 
  ggplot() + 
  geom_line(aes(x = age, y = angle, group = interaction(roi1, roi2, sex), colour = sex),
            size = .1) +
  facet_wrap(~ rsn1) + 
  theme_bw()

pairwise.angles.df %>% 
  dplyr::mutate(rsn1 == rsn2) %>% 
  ggplot() + 
  geom_point(aes(x = age, y = angle, colour = sex),
             size = .1) +
  facet_wrap(~ rsn1) + 
  theme_bw()

train.ind <- sample(c(TRUE, FALSE), nrow(ase.stats.df), replace = TRUE)
train.df <- ase.stats.df[train.ind, ]
test.df <- ase.stats.df[!train.ind, ]
y.test <- test.df$age

lin.mod <- lm(age ~ 1, data = train.df) %>% 
  step(direction = 'both',
       scope = ~ a.l + b.l + a.r + b.r + angle.diff + radius.l + radius.r + sex)
# lin.mod <- lm(age ~ 1, data = train.df) %>% 
#   step(direction = 'both',
#        scope = ~ deg.within + deg.between + sex)
summary(lin.mod)
y.hat <- predict(lin.mod, newdata = test.df)
plot(y.hat, y.test, asp = 1); abline(0, 1, col = 2)
sqrt(mean((y.hat - y.test) ^ 2))
1 - var(y.hat - y.test) / var(y.test)

inputs <- c('a.l', 'b.l', 'a.r', 'b.r', 
            'angle.diff', 'radius.l', 'radius.r')
inputs <- c('deg.within', 'deg.between')
# inputs <- c('angle.diff')
y.var <- 'age'

x.train <- train.df %>% 
  dplyr::select(inputs) %>% 
  as.matrix() %>% 
  scale() %>% 
  cbind(sex = as.numeric(as.factor(train.df$sex)) - 1)
x.test <- test.df %>% 
  dplyr::select(inputs) %>% 
  as.matrix() %>% 
  scale() %>% 
  cbind(sex = as.numeric(as.factor(test.df$sex)) - 1)
y.train <- train.df[[y.var]]
y.test <- test.df[[y.var]]
dtrain <- xgboost::xgb.DMatrix(x.train, label = y.train)
dtest <- xgboost::xgb.DMatrix(x.test, label = y.test)
xgb.mod <- xgboost::xgb.train(early_stopping_rounds = 1000, 
                              # maximize = TRUE,
                              data = dtrain,
                              nrounds = 1000,
                              params = list(nthread = 10,
                                            # objective = 'binary:logistic',
                                            # eval_metric = 'auc',
                                            # booster = 'gblinear',
                                            max_depth = 1,
                                            subsample = .5,
                                            alpha = 0,
                                            eta = .1),
                              watchlist = list(train = dtrain,
                                               eval = dtest))
xgb.mod$evaluation_log %>% 
  dplyr::mutate(train_r2 = 1 - train_rmse ^ 2 / var(train.df[[y.var]], na.rm = TRUE),
                test_r2 = 1 - eval_rmse ^ 2 / var(test.df[[y.var]], na.rm = TRUE)) %>% 
  dplyr::filter(train_r2 > 0) %>% 
  ggplot() + 
  # geom_hline(aes(yintercept = .325)) + 
  # scale_y_log10() + 
  geom_line(aes(x = iter, y = train_r2)) + 
  geom_line(aes(x = iter, y = test_r2), colour = 'red')

xgboost::xgb.importance(model = xgb.mod) %>% 
  xgboost::xgb.plot.importance(top_n = 10)

y.hat <- predict(xgb.mod, newdata = dtest)
plot(y.hat, y.test, asp = 1); abline(0, 1, col = 2)
1 - var(y.hat - y.test) / var(y.test)

A.list <- lapply(ids, function(id) {
  A <- dat.NumFibers[[id]] %>% 
    as.matrix()
  A <- (A + t(A)) / 2
  diag(A) <- 0
  A <- A / mean(A[upper.tri(A)])
  
  return(A)
})
