prognosis <- function(r) {
  df.master[df.master$rating==r,]$PD
}
loadBacktestData <- function(id){
  if (id != "all") df.defaults <- df.defaults[df.defaults$ID==id,]
  backtest <- data.frame("bank" = id, "rating"=levels(as.factor(df.defaults$Rating)))
  backtest$goods <- sapply(levels(as.factor(df.defaults$Rating)), function(x) nrow(df.defaults[df.defaults$Rating==x&df.defaults$Default==0,]))
  backtest$bads <- sapply(levels(as.factor(df.defaults$Rating)), function(x) nrow(df.defaults[df.defaults$Rating==x&df.defaults$Default==1,]))
  backtest$total <- sapply(levels(as.factor(df.defaults$Rating)), function(x) nrow(df.defaults[df.defaults$Rating==x,]))
  checkDataSize <- function(i){
    r<-as.character(backtest$rating)[i]
    Nmin <- subset(df.master, rating==r)$Nmin
    N <- backtest$total[i]
    if (N > Nmin) "N"
    else "B"
  }
  backtest <- data.frame(backtest, "model"=sapply(1:nrow(backtest), checkDataSize))
  backtest
}

# Backtesting assuming binomial distribution 
binomialTest <- function(id){
  r <- loadBacktestData(id)
  Ngoods <- sum(r$goods)
  Nbads <- sum(r$bads)
  r <- cbind(r,"PD"=sapply(as.character(r$rating),prognosis))
  r <- cbind(r,"DR"=r$bads/r$total)
  r$DR[is.na(r$DR)]<-0
  r <- cbind(r,"B"=sapply(1:nrow(r),function(i) binomialBound(r$PD[i],r$bads[i],r$total[i])))
  r$B[is.na(r$B)]<- 0
  r <- data.frame(r,"result_binom_u"=rep("G",nrow(r)))
  r <- data.frame(r,"result_binom_o"=rep("G",nrow(r)))
  r$result_binom_u<- as.character(r$result_binom_u)
  r$result_binom_o<- as.character(r$result_binom_o)
  for (i in 1:nrow(r)) {
    if (as.numeric(r$B[i]) > 0.95) r$result_binom_u[i] <- "Y"
    if (as.numeric(r$B[i]) < 0.05) r$result_binom_o[i] <- "Y"
    if (as.numeric(r$B[i]) > 0.999) r$result_binom_u[i] <- "R"
    if (as.numeric(r$B[i]) < 0.001) r$result_binom_o[i] <- "R"
  }
  r
}

# Backtesting assuming normal distribution 
normalTest <- function(id){
  r <- loadBacktestData(id)
  Ngoods <- sum(r$goods)
  Nbads <- sum(r$bads)
  r <- cbind(r,"PD"=sapply(as.character(r$rating),prognosis))
  r <- cbind(r,"DR"=r$bads/r$total)
  r$DR[is.na(r$DR)]<-0
  r$DR[is.na(r$DR)]<-0
  r <- cbind(r, "A"=r$DR-r$PD)
  r <- cbind(r,"B0.95"=sapply(1:nrow(r),function(i) normalBound(0.95,r$PD[i],r$total[i])))
  r <- cbind(r,"B0.999"=sapply(1:nrow(r),function(i) normalBound(0.999,r$PD[i],r$total[i])))
  r <- data.frame(r,"result_norm_u"=rep("G",nrow(r)))
  r$result_norm_u <- as.character(r$result_norm_u)
  r <- data.frame(r,"result_norm_o"=rep("G",nrow(r)))
  r$result_norm_o <- as.character(r$result_norm_o)
  for (i in 1:nrow(r)) {
    if (r$A[i] > r$B0.95[i]) r$result_norm_u[i] <- "Y"
    if (r$A[i] > r$B0.999[i]) r$result_norm_u[i] <- "R"
    if (-r$A[i] > r$B0.95[i]) r$result_norm_o[i] <- "Y"
    if (-r$A[i] > r$B0.999[i]) r$result_norm_o[i] <- "R"
  }
  r
}

# Backtesting assuming normal distribution + correlation
doCorrelationBacktest <- function(id){
maxDefs <- function(r,a) {
  p   <- prognosis(r)
  n   <- subset(backtest, rating == r)$total
  critValue(p,n,a,rho)
}

critValue <- function(p, n, a, rho) {
  xa   <- qnorm(1-a)
  c    <- qnorm(p)
  frac <- (c - sqrt(rho) * xa) / sqrt(1 - rho)
  Bphi <- pnorm(frac)
  phi  <- dnorm(frac)
  ceiling(n * pnorm(frac) + 0.5 * (2 * Bphi - 1 - Bphi/phi * (1 - Bphi) * sqrt((1-rho)/rho)*xa+frac))
}

trafficLight <- function(i) {
  bads <- backtest$bads[i]
  l1 <- backtest$level1[i]
  l2 <- backtest$level2[i]
  T <- "G"
  if (bads > l1) T <- "Y"
  if (bads > l2) T <- "R"
  T
}

backtest <- loadBacktestData(id)

m1 <- as.numeric(sapply(as.character(backtest$rating), maxDefs, a=a1))
m1[is.na(m1)] <- 0
m1[m1<0] <- 0
backtest <- data.frame(backtest, "level1"=m1)
m2 <- as.numeric(sapply(as.character(backtest$rating), maxDefs, a=a2))
m2[is.na(m2)] <- 0
m2[m2<0] <- 0
backtest <- data.frame(backtest, "level2"=m2)

backtest <- data.frame(backtest,"result_corr"=sapply(1:nrow(backtest), trafficLight))
backtest
}

# Calculate AUC and Gini
ratingSep <- function(id){
  q95  <- 1.36
  q999 <- 1.63
  
  r <- loadBacktestData(id)
  Ngoods <- sum(r$goods)
  Nbads <- sum(r$bads)
  dgoods <- r$goods/Ngoods
  dbads  <- r$bads/Nbads
  cdgoods <- dgoods
  cdbads <- dbads
  for (i in (length(dgoods)-1):1) cdgoods[i] <- cdgoods[i]+cdgoods[i+1]
  for (i in (length(dbads)-1):1) cdbads[i] <- cdbads[i]+cdbads[i+1]
  
  shiftCDF <- function(p, q, cdf, N){
    c <- cdf + p * q / sqrt(N)
    c[c > 1] <- 1
    c[c < 0] <- 0
    c[1]     <- 1
    c
  }
  
  AUC <- function(cdgoods, cdbads){
    x <- cdgoods[length(cdgoods)] * cdbads[length(cdbads)] / 2
    for (i in (length(cdgoods)-1):1) {
      x <- x + (cdgoods[i]-cdgoods[i+1])*(cdbads[i]+cdbads[i+1])/2
    }
    x
  }
  
  A <- AUC(cdgoods, cdbads)
  Gini <- A * 2 -1
  A95m <- AUC(shiftCDF(+1,q95,cdgoods,Ngoods),shiftCDF(-1,q95,cdbads,Nbads))
  A95p <- AUC(shiftCDF(-1,q95,cdgoods,Ngoods),shiftCDF(+1,q95,cdbads,Nbads))
  A999m <- AUC(shiftCDF(+1,q999,cdgoods,Ngoods),shiftCDF(-1,q999,cdbads,Nbads))
  A999p <- AUC(shiftCDF(-1,q999,cdgoods,Ngoods),shiftCDF(+1,q999,cdbads,Nbads))
  data.frame("Gini"=Gini,"AUC"=A,"A0.95m"=A95m,"A0.95p"=A95p,"A0.999m"=A999m,"A0.999p"=A999p,"N"=(Ngoods+Nbads))
}

# Calculate Brier score
brier <- function(id){
  r <- loadBacktestData(id)
  Ngoods <- sum(r$goods)
  Nbads <- sum(r$bads)
  r <- cbind(r,"PD"=sapply(as.character(r$rating),prognosis))
  r <- cbind(r,"DR"=r$bads/r$total)
  r$DR[is.na(r$DR)]<-0
  r <- cbind(r,"BS"=sapply(1:nrow(r), function(i) r$DR[i] * (1-r$PD[i])^2 + (1-r$DR[i]) * r$PD[i]^2))
  as.numeric(r$total %*% r$BS / (Nbads+Ngoods))
}

# Output table with Brier score, Gini, AUC + conf. intervals
analyseRating <- function(id){
  cbind("BrierScore" = brier(id), ratingSep(id))
}

# Calculate bounds for normal and binomial test
normalBound <- function(q, p, N) {
  qnorm((q+1)/2) * sqrt(p*(1-p)/N)
}
binomialBound <- function(p, Nbad, N){
  if (Nbad == 0) return(0)
  sum <- 0
  for (i in 0:Nbad) {
    sum <- sum + choose(N,i) * p^i * (1-p)^(N-i)
  }
  sum
}
