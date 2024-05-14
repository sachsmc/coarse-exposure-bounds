library(here) 
library(causaloptim)
library(data.table)
library(ggplot2)


data_to_params <- function(data, mapping, params) {
  
  valsraw <- strsplit(gsub("p", "", params), "_")
  vals <- do.call(rbind, lapply(valsraw, \(pp) {
    
    as.numeric(c(strsplit(pp[1], "")[[1]], 
    strsplit(pp[2], "")[[1]]))
    
  }))
  colnames(vals) <- c(attr(params, "rightvars"), 
                      attr(params, "condvars"))
  
  paramout <- vector(mode = "list", length = length(params))
  for(i in 1:nrow(vals)) {
    
    dsub <- as.data.frame(data[data[[mapping[attr(params, "condvars")]]] == 
                                 vals[i , attr(params, "condvars")],])
    paramout[i]<- mean(apply(as.matrix(dsub[, mapping[attr(params, "rightvars")]]), 
          MAR = 1, FUN = \(x) all(x == vals[i, attr(params, "rightvars")])))
    
    
  }
  
  names(paramout) <- params
  paramout
  
}

bnd_funcs <- readRDS(here("data/xlev_bounds.rds"))
bnd_funcs_dm <- readRDS(here("data/xlev_bounds_continuous.rds"))

data <- readRDS(here("data/analysis_data.rds"))
bnd_risk_contam <- readRDS(here("data/twolev-contam-risk.rds"))

data[, ":="(trtn = trt == "Peanut Consumption", 
            xgrp = c(0,2,1)[as.numeric(consume_group)], 
            xgrpbin = as.numeric(X2),
            gramsbin = as.numeric(consume_group != "less than .2g")
            )]

mapping1 <- c("Z" = "trtn", "X" = "xgrp", "Y" = "outcome")
mappingbin <- c("Z"= "trtn", "X"= "xgrpbin", "Y" = "outcome")

params1 <- data_to_params(data, mapping1, bnd_funcs$b3$`p{Y(X = 1) = 1} - p{Y(X = 0) = 1}`$params)
paramsbin <- data_to_params(data, mappingbin, bnd_funcs$b2$`p{Y(X = 0) = 1}`$params)

do.call(bnd_funcs$b3$`p{Y(X = 1) = 1} - p{Y(X = 0) = 1}`$bounds_func, params1)
do.call(bnd_funcs_dm$b3ms[[3]]$bounds_func, params1[bnd_funcs_dm$b3ms[[3]]$params[1:8]])
do.call(bnd_funcs$b2$`p{Y(X = 0) = 1}`$bounds_func, paramsbin)


do.call(bnd_risk_contam[[2]], paramsbin)


### bootstrap and then m-n bootstrap

bootbounds <- NULL
for(i in 1:1000) {
  data.s <- data[sample(1:nrow(data), nrow(data), replace = TRUE)]
  params1 <- data_to_params(data.s, mapping1, bnd_funcs$b3$`p{Y(X = 1) = 1} - p{Y(X = 0) = 1}`$params)
  paramsbin <- data_to_params(data.s, mappingbin, bnd_funcs$b2$`p{Y(X = 0) = 1}`$params)
  
  res <- rbind(do.call(bnd_funcs$b3$`p{Y(X = 1) = 1} - p{Y(X = 0) = 1}`$bounds_func, params1), 
               do.call(bnd_funcs_dm$b3ms[[3]]$bounds_func, params1[bnd_funcs_dm$b3ms[[3]]$params[1:8]]),
               do.call(bnd_risk_contam[[2]], paramsbin))
  res$i <- i
  res$what <- c("welldef", "illdef", "contam")
  bootbounds <- rbind(bootbounds, 
                      res)
  
}

bootbounds <- data.table(bootbounds)

bootbounds[, .(lower = quantile(lower, .05), upper = quantile(upper, .95)), by = .(what)]


n <- nrow(data)

mj <- sort(unique(floor(c(2 * sqrt(n), 0.85 ^ (0:20) * n))))
rot <- which(mj == floor(2 * sqrt(n)))

xbarstar <- array(NA, dim = c(2, B, length(mj)))

for(j in 1:length(mj)) {
  for(i in 1:B){
    
    data.s <- data[sample(1:nrow(data), 
                         mj[j], replace = (j == length(mj))),]
    
    paramsbin <- data_to_params(data.s, mappingbin, bnd_funcs$b2$`p{Y(X = 0) = 1}`$params)
    
    restab <- 
      do.call(bnd_funcs$b2$`p{Y(X = 0) = 1}`$bounds_func, paramsbin)
    
    xbarstar[,i,j] <- unlist(restab)
    
    
  }}


mhatmat <- matrix(NA, 2)

for(j in 1:2) {
    
    mhatmat[j,] <- which.min(sapply(1:(length(mj) - 1), \(xx) {
      
      ks.test(xbarstar[j,,xx], xbarstar[j,,xx+1])$statistic
      
    }))
    
}


cilim <- matrix(NA, ncol = 2)
for(j in 1:2) {
    
    cilim[j] <- quantile(xbarstar[j, , mhatmat[j,]], if(j == 1) .05 else .95)
    
}


### estimator assuming its identifiable


pdat <- subset(data, xgrp != 2) |> as.data.frame()

library(ivtools)

## ols estimate

olsest <- lm(outcome ~ xgrp +  age + sex + race, data = pdat)
confint(olsest)
coefficients(olsest)


## now g-computation


fitY <- glm(outcome ~ xgrp +  age + sex + race, data = pdat, family = binomial)
ivglm("g", "xgrp", "outcome", 
      fitZ.L = glm(trtn ~ 1, data = pdat, family = binomial),
      fitY.LZX = fitY,
      link = "identity", data = pdat) -> ivres 

summary(ivres)
confint(ivres)

ivglm("g", "xgrp", "outcome", 
      fitZ.L = glm(trtn ~ 1, data = pdat, family = binomial),
      fitY.LZX = fitY,
      link = "logit", data = pdat) -> ivreslogit


mean(plogis(predict(fitY, type = "link") - ivreslogit$est * (pdat$xgrp - 0)))
mean(plogis(predict(fitY, type = "link") - ivreslogit$est * (pdat$xgrp - 1)))


fitYbin <- glm(outcome ~ xgrpbin +  age + sex + race, data = data, family = binomial)

ivglm("g", "xgrpbin", "outcome", 
      fitZ.L = glm(trtn ~ 1, data = as.data.frame(data), family = binomial),
      fitY.LZX = fitYbin,
      link = "logit", data = as.data.frame(data)) -> ivreslogit
mean(plogis(predict(fitYbin, type = "link") - ivreslogit$est * (data$xgrp - 0)))

bsfits <- rep(NA, 1000)
for(i in 1:length(bsfits)) {
  datas <- data[sample(1:nrow(data), nrow(data), replace = TRUE),] |> as.data.frame()
  fitYbin <- glm(outcome ~ xgrpbin +  age + sex + race, data = datas, family = binomial)

  ivglm("g", "xgrpbin", "outcome", 
      fitZ.L = glm(trtn ~ 1, data = as.data.frame(datas), family = binomial),
      fitY.LZX = fitYbin,
      link = "logit", data = as.data.frame(datas)) -> ivreslogit
  bsfits[i] <- mean(plogis(predict(fitYbin, type = "link") - ivreslogit$est * (datas$xgrp - 0)))
}
