ecoli<- read.csv("ecoli-full.dat")#target=Class
levels(ecoli$Class) <- c("no", "yes")
table(ecoli$Class)
ecoli$Chg=NULL
cleveland<- read.csv("cleveland-full.dat")#target=Class
levels(cleveland$Class) <- c("no", "yes")
table(cleveland$Class)
hepatitis<- read.csv("hepatitis-full.dat")#target=Class
table(hepatitis$Class)
hepatitis$Class= recode(hepatitis$Class,"'1'=2;'2'=1")
hepatitis$Class<-as.factor(hepatitis$Class)
levels(hepatitis$Class) <- c("no", "yes")

newthyroid<- read.csv("newthyroid-full.dat")#target=Class
table(newthyroid$Class)
levels(newthyroid$Class) <- c("no", "yes")

table(Caravan$Purchase)
levels(Caravan$Purchase) <- c("no", "yes")
hepatitis$AlkPhosphate<- NULL
hepatitis$Sgot<- NULL
hepatitis$AlkPhosphate<- NULL

data=ecoli


set.seed(1)
index = round(nrow(data)*0.2,digits=0)
test.indices = sample(1:nrow(data), index)
train=data[-test.indices,]
test=data[test.indices,]

tsmote<-function (form, data, perc.over = 2, k = 5, perc.under = 2) 
{
  tgt <- which(names(data) == as.character(form[[2]]))
  minCl <- names(which.min(table(data[[tgt]])))
  minExs <- which(data[[tgt]] == minCl)
  if (tgt < ncol(data)) {
    orig.order <- colnames(data)
    cols <- 1:ncol(data)
    cols[c(tgt, ncol(data))] <- cols[c(ncol(data), tgt)]
    data <- data[, cols]
  }
  newExs <- smote.exs(data[minExs, ], ncol(data), perc.over, 
                      k)
  if (tgt < ncol(data)) {
    newExs <- newExs[, cols]
    data <- data[, cols]
  }
  selMaj <- (1:NROW(data))[-minExs]
  newdataset <- rbind(data[minExs, ],data[selMaj, ],  newExs)
  if (tgt < ncol(data)) 
    newdataset <- newdataset[, orig.order]
  
  my_list <- list(data=newdataset, min=data[minExs, ],maj=data[selMaj, ],generate=newExs)
  return(my_list)
}

smote.exs<-function (data, tgt, N, k) 
{
  nomatr <- c()
  T <- matrix(nrow = dim(data)[1], ncol = dim(data)[2] - 1)
  for (col in seq.int(dim(T)[2])) if (class(data[[col]]) %in% 
                                      c("factor", "character")) {
    T[, col] <- as.integer(data[[col]])
    nomatr <- c(nomatr, col)
  }
  else T[, col] <- data[[col]]
  if (N < 1) {
    nT <- NROW(T)
    idx <- sample(1:nT, as.integer(N * nT))
    T <- T[idx, ]
    N <- 1
  }
  p <- dim(T)[2]
  nT <- dim(T)[1]
  ranges <- apply(T, 2, max) - apply(T, 2, min)
  nexs <- as.integer(N)
  new <- matrix(nrow = nexs * nT, ncol = p)
  for (i in 1:nT) {
    xd <- scale(T, T[i, ], ranges)
    for (a in nomatr) xd[[a]] <- xd[[a]] == 0
    dd <- drop(xd^2 %*% rep(1, ncol(xd)))
    kNNs <- order(dd)[2:(k + 1)]
    for (n in 1:nexs) {
      neig <- sample(1:k, 1)
      ex <- vector(length = ncol(T))
      difs <- T[kNNs[neig], ] - T[i, ]
      new[(i - 1) * nexs + n, ] <- T[i, ] + runif(1) * 
        difs
      for (a in nomatr) new[(i - 1) * nexs + n, a] <- c(T[kNNs[neig], 
                                                          a], T[i, a])[1 + round(runif(1), 0)]
    }
  }
  newCases <- data.frame(new)
  for (a in nomatr) newCases[[a]] <- factor(newCases[[a]], 
                                            levels = 1:nlevels(data[[a]]), labels = levels(data[[a]]))
  newCases[[tgt]] <- factor(rep(data[[1, tgt]], nrow(newCases)), 
                            levels = levels(data[[tgt]]))
  colnames(newCases) <- colnames(data)
  newCases
}

##SP
##smote
test.smote <- tsmote(Class ~ ., train, perc.over = 3)
test.smote.data=test.smote$data
test.smote.maj=test.smote$maj
test.smote.min=test.smote$min
test.smote.gen=test.smote$generate




##propensity
test.smote.data$dummy=0
test.smote.data$dummy[(nrow(train)+1):(nrow(test.smote.data))]=1

propensity=rbind(test.smote.data[test.smote.data$Class=="no",],test.smote.data[test.smote.data$dummy==1,])
test.propensity=propensity %>% select(-Class)

Match.smote=matchit(dummy ~ . ,test.propensity)
promodel=cbind(test.propensity,Match.smote$distance)
matchdata=match.data(Match.smote,drop.unmatched = FALSE)
matchdata$distance=NULL
matchdata$subclass=NULL
test.smote.min$dummy=2
test.smote.min$weights=2
prin.sp=rbind(matchdata[,1:6],test.smote.min[,1:6])
sp.out=princomp(prin.sp, cor = FALSE, scores = TRUE)



matchdata$Class=propensity$Class
sp=rbind(matchdata,test.smote.min)
sp$com1=sp.out$scores[,1]
sp$com2=sp.out$scores[,2]

plot.maj=sp[sp$dummy=="0",]
plot.min=sp[sp$dummy=="2",]
plot.smote=sp[sp$dummy=="1",]
plot.propensity=sp[sp$dummy=="0" & sp$weights=="1",]
plot.sp.maj=sp[sp$dummy=="0" & sp$weights=="0",]
plot.sp.min=rbind(sp[sp$dummy=="1",],sp[sp$dummy=="2",])

##Plot
plot(plot.maj$com1,plot.maj$com2,col=1,pch=19,xlim=c(-80,80),main="SP")
legend("bottomright",
       c("maj","min","syn","match"),
       fill=c(1,2,4,5)
)
points(plot.min$com1,plot.min$com2,col=2,pch=19)
points(plot.smote$com1,plot.smote$com2,col=4)
points(plot.propensity$com1,plot.propensity$com2,col=5,pch=19)

plot(plot.sp.maj$com1,plot.sp.maj$com2,col=1,pch=19,xlim=c(-80,80),main="SP")
points(plot.min$com1,plot.min$com2,col=2,pch=19)
points(plot.smote$com1,plot.smote$com2,col=4)
legend("bottomright",
       c("maj","min","syn"),
       fill=c(1,2,4)
)
points(plot.sp.min$com1,plot.sp.min$com2,col=2,pch=19)
legend("bottomright",
       c("maj","min"),
       fill=c(1,2)
)


##PS
pro=matchit(Class ~ .,train)
plot(pro, type = 'hist', interactive = FALSE)
promodel=cbind(train,pro$distance)
matchdata=match.data(pro)


pro.maj=matchdata[matchdata$Class=="no",]

match1=pro.maj
match1$distance=NULL
match1$subclass=NULL
match1$dummy=3
match1$weights=3

pro.maj=subset(pro.maj, select = -c(weights,distance,subclass))


x1 <- apply(train, 1, paste0, collapse = ';')
x2 <- apply(pro.maj, 1, paste0, collapse = ';')
data.ps=train[!(names(x1) %in% names(x2)),]
m.smote=table(data.ps$y)[1]
n.smote=table(data.ps$y)[2]
under.smote=m.smote/(3*n.smote)-1
smote <- tsmote(Class ~ ., data.ps, perc.over = 3,perc.under=under.smote)
train=smote$data

##PSP
##1p
pro=matchit(Class ~ .,train)
matchdata=match.data(pro)
pro.maj=matchdata[matchdata$Class=="no",]


match1=pro.maj
match1$distance=NULL
match1$subclass=NULL
match1$dummy=3
match1$weights=3

pro.maj=subset(pro.maj, select = -c(weights,distance,subclass))


x1 <- apply(train, 1, paste0, collapse = ';')
x2 <- apply(pro.maj, 1, paste0, collapse = ';')

data.ps=train[!(names(x1) %in% names(x2)),]

##2s
psp <- tsmote(Class ~ ., data.ps, perc.over = 3)
data.psp=psp$data
maj.psp=psp$maj
min.psp=psp$min


##3p
data.psp$dummy=0
data.psp$dummy[(nrow(data.ps)+1):(nrow(data.psp))]=1

propensity=rbind(data.psp[data.psp$Class=="no",],data.psp[data.psp$dummy==1,])
test.propensity=propensity %>% select(-Class)

Match.smote=matchit(dummy ~ . ,test.propensity)
promodel=cbind(test.propensity,Match.smote$distance)
matchdata=match.data(Match.smote,drop.unmatched = FALSE)
matchdata$distance=NULL
matchdata$subclass=NULL
min.psp$dummy=2
min.psp$weights=2
prin.psp=rbind(matchdata[,1:6],min.psp[,1:6],match1[,1:6])
psp.out=princomp(prin.psp, cor = FALSE, scores = TRUE)

matchdata$Class=propensity$Class
psp=rbind(matchdata,min.psp,match1)
psp$com1=psp.out$scores[,1]
psp$com2=psp.out$scores[,2]

plot.maj=psp[psp$dummy=="0" | psp$dummy=="3",]
plot.min=psp[psp$dummy=="2",]
plot.ps=psp[psp$dummy=="3",]
plot.psp.fin=psp[psp$dummy=="0",]
plot.smote=psp[psp$dummy=="1",]
plot.propensity=psp[psp$dummy=="0" & psp$weights=="1",]
plot.psp.maj=psp[psp$dummy=="0" & psp$weights=="0",]
plot.psp.min=rbind(psp[psp$dummy=="1",],sp[psp$dummy=="2",])

##Plot
plot(plot.maj$com1,plot.maj$com2,col=1,pch=19,xlim=c(-80,80),main="PS")
legend("bottomright",
       c("maj","min","match"),
       fill=c(1,2,5)
)
points(plot.min$com1,plot.min$com2,col=2,pch=19)
points(plot.ps$com1,plot.ps$com2,col=5,pch=19)

plot(plot.psp.fin$com1,plot.psp.fin$com2,col=1,pch=19,xlim=c(-80,80),main="PSP")
legend("bottomright",
       c("maj","min","syn","match"),
       fill=c(1,2,4,5)
)
points(plot.min$com1,plot.min$com2,col=2,pch=19)
points(plot.smote$com1,plot.smote$com2,col=4)
points(plot.propensity$com1,plot.propensity$com2,col=5,pch=19)


plot(plot.psp.maj$com1,plot.psp.maj$com2,col=1,pch=19,xlim=c(-80,80),main="PSP")
points(plot.min$com1,plot.min$com2,col=2,pch=19)
legend("bottomright",
       c("maj","min","syn"),
       fill=c(1,2,4)
)
points(plot.smote$com1,plot.smote$com2,col=4)
points(plot.smote$com1,plot.smote$com2,col=2,pch=19)
legend("bottomright",
       c("maj","min"),
       fill=c(1,2)
)
















matchdata=match.data(Match.smote)
pro.maj=matchdata[matchdata$dummy==0,]
pro.maj=subset(pro.maj, select = -c(weights,distance,subclass,dummy))
pro.maj$Class="no"

x1 <- apply(data.ps, 1, paste0, collapse = ';')
x2 <- apply(pro.maj, 1, paste0, collapse = ';')

train.maj=data.ps[!(names(x1) %in% names(x2)),]
train=rbind(train.maj,psp$generate)
