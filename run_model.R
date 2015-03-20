library(doParallel)
library(plyr)
library(randomForest)

nodes <- detectCores()
cl <- makeCluster(nodes)
registerDoParallel(cl)

fname = "result8_2_3_icode_new.csv.gz"
drivers = list.files("./drivers")
randomDrivers = sample(drivers, size = 20) 
write.csv(randomDrivers, paste0("randomDrivers_", fname, date() ), row.names=F)

featurext <- function(trip)
{
  euclid <- function(a) sqrt((a["x"])^2 + (a["y"])^2)
  
  dotp <- function(a,b) a["x"]*b["x"] + a["y"]*b["y"]
  
  speeds <- sqrt(diff(trip$x)^2 + diff(trip$y)^2)
  accs <- diff(speeds)
  paccs <- accs[accs>0 & accs < 10]  
  naccs <- accs[accs<0 & accs > -10]  
  speedv <- speeds[speeds<50]
  
  maxd=nrow(trip)
  if (euclid(trip[maxd,]) == 0) {maxd <- which(apply(trip,1,euclid)==max(apply(trip,1,euclid)))[1]}
  sina <- trip$x[maxd] / euclid(trip[maxd,])
  cosa <- trip$y[maxd] / euclid(trip[maxd,])
  
  nr <- nrow(trip)
  trip1 <- trip[1:(nr-2),]
  trip2 <- trip[2:(nr-1),]
  trip3 <- trip[3:nr,]
  
  angle <- acos(dotp(trip2-trip1, trip3-trip2)/(euclid(trip2-trip1)*euclid(trip3-trip2)))
  angle$x[is.nan(angle$x)]=0
  angle <- angle[angle$x < 0.5,]
  ang <- hist(angle, seq(0,0.5,0.025), plot=F)$counts/length(angle)
  
  tripr <- t(apply(as.matrix(trip), 1, function(a) {c(a["x"]*cosa - a["y"]*sina, a["x"]*sina + a["y"]*cosa)}))
  maxs <- apply(tripr,2,max)
  mins <- apply(tripr,2,min)

  finess=1
  xmax = floor(maxs[1]/finess)
  ymax = floor(maxs[2]/finess)
  xmin = floor(mins[1]/finess)
  ymin = floor(mins[2]/finess)
  
  return(c(hist(paccs, seq(0,10,0.5), plot=F)$density/2, 
           hist(naccs, seq(-10,0,0.5), plot=F)$density/2, 
           hist(speedv, seq(0,50,2), plot=F)$density,
           nrow(trip),
           sum(speeds),
           quantile(speeds, seq(0,1,0.1), names=F),
           xmax - xmin,
           ymax - ymin,
           ang,
           sum(angle)
          )
        )
}

processDriver <- function(driver){
  
  require(randomForest)
  
  dirPath = paste0("./drivers/", driver, '/')
  currentTrip <- function(i) {
    trip = read.csv(paste0(dirPath, i, ".csv"))
    return(c(featurext(trip), 1))    
  }   
  currentData = ldply(1:200, currentTrip, .parallel=T)
  colnames(currentData)[ncol(currentData)] = "t"
  
  id = which(refData[,c("driver")]==driver)
  rData = refData[,(1:ncol(refData)-1)]
  #if (length(id)>0) rData[id, c("t")]=1  
  
  train = rbind(currentData, rData)
  train = as.data.frame(train)
  train$t <- as.factor(train$t)
  g = randomForest(t ~., data = train, na.action=na.roughfix)
  p = predict(g, type = "prob", na.action=na.roughfix)
  
  data.frame(driver_trip=paste0(driver, "_", 1:200),
             prob=p[1:200,2],
             row.names=NULL,
             stringsAsFactors=F)
}

refTrips <- function(d) {
  tripData <- function(i) {
    trip = read.csv(paste0(dirPath, i, ".csv"))
    return(c(featurext(trip), 0, as.numeric(d)))    
  }
  l = length(randomDrivers)
  
  dirPath = paste0("./drivers/", d, '/')
  ixs <- sample(1:200, size=floor(800/l))
  refTripData = ldply(ixs, tripData, .parallel = T)
  colnames(refTripData)[ncol(refTripData)-1] = "t"
  data.frame(refTripData)
}

clusterExport(cl, list("featurext", "randomDrivers", "processDriver"))
refData = ldply(randomDrivers, refTrips, .parallel = T)
colnames(refData)[ncol(refData)]="driver"

clusterExport(cl, list("refData"))
submission <- ldply(drivers, processDriver, .parallel=T)

write.csv(submission,gzfile(fname), row.names=F, quote=F)
stopCluster(cl)
