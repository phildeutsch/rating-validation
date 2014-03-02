# load masterscale
filename.master <- 'master'
# load rating data
filename.data   <- 'data_input'

a1 <- 0.95
a2 <- 0.999
rho <- 0.01

# Load data
master <- read.csv(paste(filename.master,'.csv',sep=""),head=T,sep=",")
df.master<-as.data.frame(master)
names(df.master) <- c('rating','PD')
df.master <- data.frame(df.master,"Nmin"=sapply(1:nrow(df.master),function(i) ceiling(9/(df.master[i,2]*(1-df.master[i,2])))))
data <- read.csv(paste(filename.data,'.csv',sep=""),head=T,sep=",")
# Define the columns to be used
df.defaults <- data.frame(data[,1],as.character(data[,4]),data[,5],data[,33])
names(df.defaults) <- c("ID","BehaviourRating","Rating","Default")
df.defaults <- subset(df.defaults, BehaviourRating != "NR")

# Check model suitability
result <- loadBacktestData("all")
for (id in levels(as.factor(df.defaults$ID))) {
  result <- rbind(result,loadBacktestData(id))
}
result
write.table(result, paste(filename.data,'_resultModel.csv',sep=""),sep=",")

# Calculate Brier, Gini, AUC + Conf.Interval
result <- cbind("bank"="all",analyseRating("all"))
for (id in levels(as.factor(df.defaults$ID))) {
  result <- rbind(result,cbind("bank"=id,analyseRating(id)))
}
write.table(result, paste(filename.data,'_resultBGA.csv',sep=""),sep=",")

# Normal Backtest for all banks
result <- normalTest("all")
for (id in levels(as.factor(df.defaults$ID))) {
  result <- rbind(result,normalTest(id))
}
write.table(result, paste(filename.data,'_resultNormal.csv',sep=""),sep=",")

# Binomial Backtest for all banks
result <- binomialTest("all")
for (id in levels(as.factor(df.defaults$ID))) {
  result <- rbind(result,binomialTest(id))
}
write.table(result, paste(filename.data,'_resultBinomial.csv',sep=""),sep=",")

# Correlation Backtest for all banks
result <- doCorrelationBacktest("all")
for (id in levels(as.factor(df.defaults$ID))) {
  result <- rbind(result,doCorrelationBacktest(id))
}
write.table(result, paste(filename.data,'_resultCorrelation_noNR.csv',sep=""),sep=",")
