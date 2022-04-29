

algorithm=function(data,target,method="PSP",matching.method="Cutpoint",point=0.2,interval=1.645){
  set.seed(987)
  ##
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
  
  matching.ps=function(data,target="y",method,point,interval){
    #########initial setting#####################
    
    names(data)[names(data) == target] <-"dummy"
    data$dummy=as.factor(data$dummy)
    level=levels(data$dummy)
    table=table(data$dummy)
    if (table[1]>=table[2]){
      min=data[data$dummy==level[2],]
      maj=data[data$dummy==level[1],]
    } else {
      min=data[data$dummy==level[1],]
      maj=data[data$dummy==level[2],]
    }
    n=nrow(min)
    list=list()
    
    ############Nearest###########################
    
    if (method=="Nearest") {
      
      for (i in 1:n){
        
        list[i]=which(abs(maj$distance-min$distance[i])==min(abs(maj$distance-min$distance[i])))
        
      }
      unique.list=unique(list)
      list.matched=sort(unlist(unique.list))
      matched=maj[list.matched,]
      list.all=c(1:nrow(maj))
      unmatched=maj[!(list.all %in% list.matched),]
    }
    else if (method=="Cutpoint"){
      
      matched=maj[maj$distance>point,]
      unmatched=maj[maj$distance<point,]
      
    }
    else {
      
      for (i in 1:n){
        
        list[i]=which(abs(maj$distance-min$distance[i])==min(abs(maj$distance-min$distance[i])))
        std=sd(min$distance)
        t=maj$distance[unlist(list[i])]
        l=min$distance[i]-interval*std/n
        u=min$distance[i]+interval*std/n
        if ( l<t & t<u){
          maj$distance[c(unlist(list[i]))]=-1
        }  else {
          list[i]=NULL
        }
        
      }
      unique.list=unique(list)
      list.matched=sort(unlist(unique.list))
      matched=maj[list.matched,]
      list.all=c(1:nrow(maj))
      unmatched=maj[!(list.all %in% list.matched),]
      
    }
    
    
    if (target=="dummy"){
      matched$dummy=NULL
      unmatched$dummy=NULL
      min$dummy=NULL
    }  else  {
      names(matched)[names(matched) == "dummy"] <-target
      names(unmatched)[names(unmatched) == "dummy"] <-target
      names(min)[names(min) == "dummy"] <-target
    }
    
    matched$distance=NULL
    unmatched$distance=NULL  
    min$distance=NULL
    matchingdata=rbind(unmatched,min)
    matching_list <- list(matched=matched,unmatched=unmatched,data=matchingdata,list)
    return(matching_list)
    
  }
  
  
  
  names(data)[names(data) == target] <-"y"  #rename target
  
  train=data
  ##
  
  if (method=="origin"){
    ##set data to training and testing
    
    train=data
    
    
    ##
    
  }
  else if (method=="smote"){
    m.smote=table(train$y)[1]
    n.smote=table(train$y)[2]
    over.smote=m.smote/(2*n.smote)
    set.seed(1234)
    smote <- smote(y ~ ., data, perc.over = 3)
    train=smote
    
    
  }
  else if (method=="PS"){
    
    pro=matchit(y ~ .,train)
    distance=pro$distance
    psdata=cbind(train,distance)
    matchdata=matching.ps(psdata,target="y",matching.method,point,interval)
    data.ps=matchdata$data
    
    m.smote=table(train$y)[1]
    n.smote=table(train$y)[2]
    over.smote=m.smote/(2*n.smote)
    
    smote <- smote(y ~ ., data.ps, perc.over = 3)
    train=smote
    
  }
  else if (method=="SP"){
    ##smote
    
    m.smote=table(train$y)[1]
    n.smote=table(train$y)[2]
    over.smote=m.smote/(2*n.smote)
    test.smote <- tsmote(y ~ ., train, perc.over = 3)
    test.smote.data=test.smote$data
    test.smote.maj=test.smote$maj
    test.smote.min=test.smote$min
    test.smote.gen=test.smote$generate
    
    ##propensity
    test.smote.data$dummy=0
    test.smote.data$dummy[(nrow(train)+1):(nrow(test.smote.data))]=1
    
    propensity=rbind(test.smote.data[test.smote.data$y=="no",],test.smote.data[test.smote.data$dummy==1,])
    test.propensity=propensity %>% select(-y)
    
    Match.smote=matchit(dummy ~ . ,test.propensity)
    distance=Match.smote$distance
    spdata=cbind(propensity,distance)
    matchdata=matching.ps(spdata,target="dummy",matching.method,point,interval)
    SP=rbind(matchdata$data,test.smote.min)
    train.maj=SP[SP$y=="no",]
    train.min=SP[SP$y=="yes",]
    n.train.min=nrow(train.min)
    if(nrow(train.maj)>3/2*n.train.min){
      train.maj.sample=train.maj[sample(nrow(train.maj), 3/2*n.train.min),]
    } else{
      train.maj.sample=train.maj
    }
    train=rbind(train.maj.sample,train.min)
    
    
  }
  else if (method=="PSP"){
    
    ##1p
    pro=matchit(y ~ .,train)
    distance=pro$distance
    psdata=cbind(train,distance)
    matchdata=matching.ps(psdata,target="y",matching.method,point,interval)
    data.ps=matchdata$data
    
    ##2s
    m.smote=table(data.ps$y)[1]
    n.smote=table(data.ps$y)[2]
    over.smote=m.smote/(2*n.smote)
    psp <- tsmote(y ~ ., data.ps, perc.over = 3)
    data.psp=psp$data
    maj.psp=psp$maj
    min.psp=psp$min
    
    
    ##3p
    data.psp$dummy=0
    data.psp$dummy[(nrow(data.ps)+1):(nrow(data.psp))]=1
    
    propensity=rbind(data.psp[data.psp$y=="no",],data.psp[data.psp$dummy==1,])
    test.propensity=propensity %>% select(-y)
    
    Match.smote=matchit(dummy ~ . ,test.propensity)
    distance=Match.smote$distance
    spdata=cbind(propensity,distance)
    matchdata=matching.ps(spdata,target="dummy",matching.method,point,interval)
    PSP=rbind(matchdata$data,min.psp)
    train.maj=PSP[PSP$y=="no",]
    train.min=PSP[PSP$y=="yes",]
    n.train.min=nrow(train.min)
    if(nrow(train.maj)>3/2*n.train.min){
      train.maj.sample=train.maj[sample(nrow(train.maj), 3/2*n.train.min),]
    } else{
      train.maj.sample=train.maj
    }
    train=rbind(train.maj.sample,train.min)
    
    
    
  }
  else if (method=="Tomek"){
    ##Tomek##
    
    
    smote <- tsmote(y ~ ., train, perc.over = 3,perc.under=5)
    train=smote$data
    smote.x = as.matrix(train %>% select(-y))
    smote.y = as.matrix(recode(train$y,"'no'=0;'yes'=1"))
    tomek <- ubTomek(smote.x,smote.y)
    tomek$X=as.data.frame( tomek$X)
    tomek$X[] <- lapply(  tomek$X, as.numeric)
    tomek.train<-cbind( tomek$X, tomek$Y)
    colnum=length(tomek.train)
    colnames(tomek.train)[colnum] <- "y"
    levels (tomek.train$y) <- c("no", "yes")
    train<-tomek.train
    
    
    
    
    
    
  }
  else if (method=="enn"){
    ##enn##
    smote <- tsmote(y ~ ., train, perc.over = 3,perc.under=5)
    train=smote$data
    smote.x = as.matrix(train %>% select(-y))
    smote.y = as.matrix(recode(train$y,"'no'=0;'yes'=1"))
    enn <-ubENN(smote.x, smote.y, k = 3, verbose = TRUE)
    enn.train <-cbind(as.data.frame(enn$X), as.data.frame(enn$Y))
    colnum=length(enn.train)
    colnames(enn.train)[colnum] <- "y"
    levels (enn.train$y) <- c("no", "yes")
    train<-enn.train
    
    
    
  }  
  else if (method=="ROSE"){
    rose=ROSE(y~.,train)
    train= rose$data
    
    
  }
  ##Plot##
  
  match=matchit(y ~ . ,train)
  
  pscore=match$distance
  Class=train$y
  plot=plot(Class,pscore,xlab="Class",ylab="PScore",outline=FALSE,ylim=c(0,1))
  title(main = method)
  
  ##PCAPlot##
  
  
  
  names(train)[names(train) == "y"] <-target
  my_list <- list(data=train,plot=plot,match=matchdata)
  return(my_list)
  
}

Measure=function(train,test,target,ROC=TRUE,cv="kfolds",folds=5){
  names(train)[names(train) == target] <-"y"  
  names(test)[names(test) == target] <-"y"  
  
  ##Measurement##
  Measure=function(predict,y){
    tp <- length(which(ifelse(predict == "yes" & y == "yes", TRUE, FALSE) == TRUE))
    tn <- length(which(ifelse(predict == "no" & y == "no", TRUE, FALSE) == TRUE))
    fp <- length(which(ifelse(predict == "yes" & y == "no", TRUE, FALSE) == TRUE))
    fn <- length(which(ifelse(predict == "no" & y == "yes", TRUE, FALSE) == TRUE))
    tprate=tp/(tp+fn)
    tnrate=tn/(tn+fp)
    ppvalue=tp/(tp+fp)
    
    
    
    
    gmean=round(sqrt(tprate*tnrate),4)
    fmeasure=round(2*tprate*ppvalue/(tprate+ppvalue),4)
    N=tn+tp+fn+fp
    S=(tp+fn)/N
    P=(tp+fp)/N
    MCC=round((tp/N-S*P)/sqrt(P*S*(1-S)*(1-P)),4)
    CK=round(2*(tp*tn-fp*fn)/((tp+fp)*(fp+tn)+(fn+tn)*(tp+fn)),4)
    return(c("gmean=",gmean,"fmeasure=",fmeasure,"MCC=",MCC,"Cohen's kappa ",CK))
  }
  
  
  
  if(ROC==TRUE){
    ##set records##
    #    records = matrix(NA, nrow=2, ncol=5)
    #    colnames(records) <- c("gmean","fmeasure","ROC","MCC","CK") 
    #    rownames(records) <- c("Logistic","Random Forests")
    records = matrix(NA, nrow=5, ncol=5)
    colnames(records) <- c("gmean","fmeasure","ROC","MCC","CK") 
    rownames(records) <- c("Logistic","Tree","KNN","Random Forests","Neural Network")
    if (cv=="NO"){
      train.control <- trainControl(method = "none")
      # Train the model
      ##glm
      glm <- caret::train(y ~., data= train, method = "glm",
                          trControl = train.control)
      
      pred_glm.train = predict(glm,train,type="raw")
      pred_glm.test = predict(glm, test,type="raw")
      measure.glm=Measure(pred_glm.train,train$y)
      measure.glm.t=Measure(pred_glm.test,test$y)
      roc_glm.train <- roc(as.numeric(pred_glm.train),as.numeric(train$y))
      auc_glm.train <- as.numeric(auc(roc_glm.train))
      roc_glm.test <- roc(as.numeric(pred_glm.test),as.numeric(test$y))
      auc_glm.test <- as.numeric(auc(roc_glm.test))    
      records[1,] <- c(measure.glm.t[2],measure.glm.t[4],auc_glm.test,measure.glm.t[6],measure.glm.t[8])#write tree into the second row
      
      ##rf
      RF <- caret::train(y ~., data= train, method = "rf",
                         trControl = train.control)
      pred_RF.train=predict(RF,train,type="raw")
      pred_RF.test=predict(RF,test,type="raw")
      measure.rf=Measure(pred_RF.train,train$y)
      measure.rf.t=Measure(pred_RF.test,test$y)
      roc_RF.train <- roc(as.numeric(pred_RF.train),as.numeric(train$y))
      auc_RF.train <- as.numeric(auc(roc_RF.train))
      roc_RF.test <- roc(as.numeric(pred_RF.test),as.numeric(test$y))
      auc_RF.test <- as.numeric(auc(roc_RF.test))
      records[2,] <- c(measure.rf.t[2],measure.rf.t[4],auc_RF.test,measure.rf.t[6],measure.rf.t[8])#write rf into the 4th row
      
    } else {
      
      #     cvIndex <- createFolds(factor(train$y), folds, returnTrain = T)
      train.control <- trainControl(method = 'cv', number = folds,savePredictions = TRUE, 
                                    classProbs = TRUE, 
                                    verboseIter = TRUE)
      
      
      
      # Train the model
      ##glm
      glm <- caret::train(y ~., data= train, method = "glm",
                          trControl = train.control)
      
      
      pred_glm.test = predict(glm, test,type="raw")
      
      measure.glm.t=Measure(pred_glm.test,test$y)
      pred_glm.test[1]="yes"      
      roc_glm.test <- roc(as.numeric(pred_glm.test),as.numeric(test$y))
      auc_glm.test <- as.numeric(auc(roc_glm.test))    
      records[1,] <- c(measure.glm.t[2],measure.glm.t[4],auc_glm.test,measure.glm.t[6],measure.glm.t[8])#write tree into the second row
      
      ##tree
      tree <- caret::train(y ~., data= train, method = "rpart",
                           trControl = train.control)
      
      
      pred_tree.test = predict(tree, test,type="raw")
      measure.tree.t=Measure(pred_tree.test,test$y)
      pred_tree.test[1]="yes"
      roc_tree.test <- roc(as.numeric(pred_tree.test),as.numeric(test$y))
      auc_tree.test <- as.numeric(auc(roc_tree.test))   
      
      records[2,] <- c(measure.tree.t[2],measure.tree.t[4],auc_tree.test,measure.tree.t[6],measure.tree.t[8])
      
      ##KNN
      knn <- caret::train(y ~., data= train, method = "knn",
                          trControl = train.control)      
      
      pred_knn.test = predict(knn, test,type="raw")
      measure.knn.t=Measure(pred_knn.test,test$y)
      pred_knn.test[1]="yes"
      roc_knn.test <- roc(as.numeric(pred_knn.test),as.numeric(test$y))
      auc_knn.test <- as.numeric(auc(roc_knn.test))   
      
      records[3,] <- c(measure.knn.t[2],measure.knn.t[4],auc_knn.test,measure.knn.t[6],measure.knn.t[8])
      
      ##rf
      RF <- caret::train(y ~., data= train, method = "rf",
                         trControl = train.control)
      
      pred_RF.test=predict(RF,test,type="raw")
      
      measure.rf.t=Measure(pred_RF.test,test$y)
      pred_RF.test[1]="yes"
      roc_RF.test <- roc(as.numeric(pred_RF.test),as.numeric(test$y))
      auc_RF.test <- as.numeric(auc(roc_RF.test))
      records[4,] <- c(measure.rf.t[2],measure.rf.t[4],auc_RF.test,measure.rf.t[6],measure.rf.t[8])#write rf into the 4th row
      #Neural Network
      nn <- caret::train(y ~., data= train, method = "nnet",
                         trControl = train.control)      
      
      pred_nn.test = predict(nn, test,type="raw")
      measure.nn.t=Measure(pred_nn.test,test$y)
      pred_nn.test[1]="yes"
      roc_nn.test <- roc(as.numeric(pred_nn.test),as.numeric(test$y))
      auc_nn.test <- as.numeric(auc(roc_nn.test))   
      
      records[5,] <- c(measure.nn.t[2],measure.nn.t[4],auc_nn.test,measure.nn.t[6],measure.nn.t[8])
    }
    
  } else{
    ##set records##
    records = matrix(NA, nrow=5, ncol=4)
    colnames(records) <- c("train.gmean","train.fmeasure","test.gmin","test.fmeasure") 
    rownames(records) <- c("Logistic","Tree","KNN","Random Forests","SVM")
    
    if (cv=="Kfolds"){
      train.control <- trainControl(method = "cv", number = 5)
      # Train the model
      ##glm
      glm <- caret::train(y ~., data= train, method = "glm",
                          trControl = train.control)
      
      pred_glm.train = predict(glm,train,type="raw")
      pred_glm.test = predict(glm, test,type="raw")
      measure.glm=Measure(pred_glm.train,train$y)
      measure.glm.t=Measure(pred_glm.test,test$y)
      records[1,] <- c(measure.glm[2],measure.glm[4],measure.glm.t[2],measure.glm.t[4])#write tree into the second row
      
      ##tree
      tree <- caret::train(y ~., data= train, method = "rpart",
                           trControl = train.control)
      
      pred_tree.train = predict(tree,train,type="raw")
      pred_tree.test = predict(tree, test,type="raw")
      measure.tree=Measure(pred_tree.train,train$y)
      measure.tree.t=Measure(pred_tree.test,test$y)
      
      records[2,] <- c(measure.tree[2],measure.tree[4],measure.tree.t[2],measure.tree.t[4])#write tree into the second row
      ##knn
      knn <- caret::train(y ~., data= train, method = "knn",
                          trControl = train.control)
      pred_knn.train=predict(knn,train,type="raw")
      pred_knn.test=predict(knn,test,type="raw")
      measure.knn=Measure(pred_knn.train,train$y)
      measure.knn.t=Measure(pred_knn.test,test$y)
      
      records[3,] <- c(measure.knn[2],measure.knn[4],measure.knn.t[2],measure.knn.t[4])#write knn into the third row
      ##rf
      RF <- caret::train(y ~., data= train, method = "rf",
                         trControl = train.control)
      pred_RF.train=predict(RF,train,type="raw")
      pred_RF.test=predict(RF,test,type="raw")
      measure.rf=Measure(pred_RF.train,train$y)
      measure.rf.t=Measure(pred_RF.test,test$y)
      
      records[4,] <- c(measure.rf[2],measure.rf[4],measure.rf.t[2],measure.rf.t[4])#write rf into the 4th row
      
      ##svm
      svm=caret::train(y ~., data= train, method = "svmLinear",
                       trControl = train.control)
      pred_svm.train = predict(svm,train,type="raw")
      pred_svm.test = predict(svm,test,type="raw")
      measure.svm=Measure(pred_svm.train,train$y)
      measure.svm.t=Measure(pred_svm.test,test$y)
      
      records[5,] <- c(measure.svm[2],measure.svm[4],measure.svm.t[2],measure.svm.t[4])#write svm into the 5th row
    }
    else if (cv=="NO"){
      train.control <- trainControl(method = "none")
      # Train the model
      ##glm
      glm <- caret::train(y ~., data= train, method = "glm",
                          trControl = train.control)
      
      pred_glm.train = predict(glm,train,type="raw")
      pred_glm.test = predict(glm, test,type="raw")
      measure.glm=Measure(pred_glm.train,train$y)
      measure.glm.t=Measure(pred_glm.test,test$y)
      records[1,] <- c(measure.glm[2],measure.glm[4],measure.glm.t[2],measure.glm.t[4])#write tree into the second row
      
      ##tree
      tree <- caret::train(y ~., data= train, method = "J48",
                           trControl = train.control)
      
      pred_tree.train = predict(tree,train,type="raw")
      pred_tree.test = predict(tree, test,type="raw")
      measure.tree=Measure(pred_tree.train,train$y)
      measure.tree.t=Measure(pred_tree.test,test$y)
      
      records[2,] <- c(measure.tree[2],measure.tree[4],measure.tree.t[2],measure.tree.t[4])#write tree into the second row
      ##knn
      knn <- caret::train(y ~., data= train, method = "knn",
                          trControl = train.control)
      pred_knn.train=predict(knn,train,type="raw")
      pred_knn.test=predict(knn,test,type="raw")
      measure.knn=Measure(pred_knn.train,train$y)
      measure.knn.t=Measure(pred_knn.test,test$y)
      
      records[3,] <- c(measure.knn[2],measure.knn[4],measure.knn.t[2],measure.knn.t[4])#write knn into the third row
      ##rf
      RF <- caret::train(y ~., data= train, method = "rf",
                         trControl = train.control)
      pred_RF.train=predict(RF,train,type="raw")
      pred_RF.test=predict(RF,test,type="raw")
      measure.rf=Measure(pred_RF.train,train$y)
      measure.rf.t=Measure(pred_RF.test,test$y)
      
      records[4,] <- c(measure.rf[2],measure.rf[4],measure.rf.t[2],measure.rf.t[4])#write rf into the 4th row
      
      ##svm
      svm=train(y ~., data= train, method = "svmLinear",
                trControl = train.control)
      pred_svm.train = predict(svm,train,type="raw")
      pred_svm.test = predict(svm,test,type="raw")
      measure.svm=Measure(pred_svm.train,train$y)
      measure.svm.t=Measure(pred_svm.test,test$y)
      
      records[5,] <- c(measure.svm[2],measure.svm[4],measure.svm.t[2],measure.svm.t[4])#write svm into the 5th row
    }
  }
  table=table(train$y)
  
  my_list <- list(records=records,table=table)
  return(my_list)
  
}
