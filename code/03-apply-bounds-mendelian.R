library(here)
library(causaloptim)
library(data.table)
library(xtable)

rdat <- read.table(here("rawdata/mendelian-data.txt"), sep = "\t", header = TRUE)

Xgroup <- rdat[, 1]
Xmap1 <- c(0, 2, 2, 1, 1)
Xmap2 <- c(0, 2, 3, 1, 1)
names(Xmap1) <- names(Xmap2) <- Xgroup

dat <- do.call(rbind, lapply(rdat[, 2:4], \(col) {
  data.frame(Y = c(0, 1), count = as.numeric(unlist(strsplit(col, "/"))))
}))
dat$Z <- rep(c(0, 1, 2), each = 2 * 5)
dat$Xraw <- rep(Xgroup, each = 2)
dat$X1 <- Xmap1[dat$Xraw]
dat$X2 <- Xmap2[dat$Xraw]

Ztots <- by(dat$count, dat$Z, sum)
dat$Ztot <- NA
dat$Ztot[dat$Z == 0] <- Ztots[[1]]
dat$Ztot[dat$Z == 1] <- Ztots[[2]]
dat$Ztot[dat$Z == 2] <- Ztots[[3]]

datsum <- data.table(dat)

datsum[, jack := sprintf("%.0f (%.1f%%)", count, 100 * count / Ztot)]
datsum[, Xfact := factor(Xraw, levels = Xgroup, ordered = TRUE)]
dcast(datsum[Y == 0], Xfact ~ Z, value.var = "jack") |> xtable() |> print(include.rownames = FALSE)
dcast(datsum[Y == 1], Xfact ~ Z, value.var = "jack") |> xtable() |> print(include.rownames = FALSE)

prob_1 <- by(dat[, c("count", "Ztot")], list(Y = dat$Y, X = dat$X1, Z = dat$Z), 
             \(df) {
               data.frame(count = sum(df$count), prob = sum(df$count) / df$Ztot[1])
             }) |> array2DF()
prob_2 <- by(dat[, c("count", "Ztot")], list(Y = dat$Y, X = dat$X2, Z = dat$Z), 
                 \(df) {
                   data.frame(count = sum(df$count), prob = sum(df$count) / df$Ztot[1])
                 }) |> array2DF()

prob_1$names <- with(prob_1, sprintf("p%s%s_%s", X, Y, Z))
prob_2$names <- with(prob_2, sprintf("p%s%s_%s", X, Y, Z))

### running the bounds

trebnds <- readRDS(here("data/threeiv-threeX-classic.rds"))$boundsFunction
treprobs <- sapply(names(formals(trebnds)), \(probN) {
  res <- subset(prob_1, names == probN)$prob
  res
}, simplify = FALSE)

do.call(trebnds, treprobs)


tre.illbnds <- readRDS(here("data/threeiv-threeX-contam.rds"))[[2]]
do.call(tre.illbnds, treprobs)


tre.contbnds <- readRDS(here("data/threeiv-threeX-illdef.rds"))[[2]]

treprobs2 <- sapply(names(formals(tre.contbnds)), \(probN) {
  res <- subset(prob_1, names == probN)$prob
  res
}, simplify = FALSE)
do.call(tre.contbnds, treprobs2)

## parametric bootstrap

bootprobs <- function(prob) {
  
  prob <- data.table(prob)
  prob[, newprob := rmultinom(1, sum(count), prob) / sum(count), by = .(Z)]
 
  prob
  
}


boottre <- NULL
for(i in 1:1000) {
  
  prob_s <- bootprobs(prob_1)
  plist_s <- sapply(names(formals(trebnds)), \(probN) {
    res <- subset(prob_s, names == probN)$newprob
    res
  }, simplify = FALSE)
  
  boottre <- rbind(boottre, do.call(trebnds, plist_s))
  
  
}

quantile(boottre$lower, .05)
quantile(boottre$upper, .95)



quadbnds <- readRDS(here("data/threeiv-fourX-classic.rds"))$boundsFunction
quadprobs <- sapply(names(formals(quadbnds)), \(probN) {
  res <- subset(prob_2, names == probN)$Value
  res
}, simplify = FALSE)

do.call(quadbnds, quadprobs)


quad.illbnds <- readRDS(here("data/threeiv-fourX-contam.rds"))[[2]]
do.call(quad.illbnds, quadprobs)


quad.contbnds <- readRDS(here("data/threeiv-fourX-illdef.rds"))[[2]]

quadprobs2 <- sapply(names(formals(quad.contbnds)), \(probN) {
  res <- subset(prob_2, names == probN)$Value
  res
}, simplify = FALSE)
do.call(quad.contbnds, quadprobs2)



## simple wald method

datsum2 <- datsum[X2 %in% c(0,1)]

dplu <- datsum2[, .(count = sum(count)), by = .(Z, X1, Y)]
dplu[, probX := rep(sum(count[X1 == 1]) / sum(count), .N), by = .(Z)]

ndplu <- dplu[, .(Yn = rep(Y, count), X = rep(probX, count)), by = .(Z, X1, Y)]

lm(Yn ~ X, data = ndplu)
lm(Yn ~ X1, data = ndplu) |> confint()  ## ols

## now bootstrap this

dplu[, mprob := count / sum(count), by = .(Z)]

bres <- rep(NA, 1000)
for(i in 1:length(bres)){
  dplus <- dplu[, countstar := c(rmultinom(1, sum(count), mprob)), by = .(Z)]
  dplus[, probX := rep(sum(count[X1 == 1]) / sum(countstar), .N), by = .(Z)]
  
  ndplu <- dplus[, .(Yn = rep(Y, countstar), X = rep(probX, countstar)), by = .(Z, X1, Y)]
  
  bres[i] <- lm(Yn ~ X, data = ndplu)$coefficients[2]
  
  
}

quantile(bres, c(0.025, 0.975))
