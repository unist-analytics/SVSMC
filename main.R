libs<-c("tm","plyr","class","kknn","UBL","mclust","zeallot","FSelector","GA","caret","ModelMetrics","mltools")
lapply(libs,require,character.only=TRUE)
options(stringsAsFactors = FALSE)
#path<- "../"
#setwd(paste(path))


source('source.R', echo=FALSE)

load(file="tdm.RData")



## start 

beta_bound <- c()
beta_bound3 <- c()

maucs<-c()
Gmetrics<-c()
dans<-c()


for (nn in 1:2){
  mauc_s<-data.frame(matrix(NA,nrow=1,ncol=8))
  names(mauc_s)<-c("ori","syn","syn3","smo","scut","edit","IncuBAtE")
  Gmetric_s<-data.frame(matrix(NA,nrow=1,ncol=8))
  names(Gmetric_s)<-c("ori","syn","syn3","smo","scut","edit","IncuBAtE")
  dan_s<-data.frame(matrix(NA,nrow=1,ncol=8))
  names(dan_s)<-c("ori","syn","syn3","smo","scut","edit","IncuBAtE")
  
  
  c(train.kknn,test.kknn) %<-% c(0,0)
  
  tdm = c()
  for (i in 1:5){
    Tdm =c()
    Tdm <- list(subset(tdm.stack,targetCandidate==i))
    tdm <- c(tdm,Tdm)
  }
  
  for (i in 1:5){
    train.idx = c()
    test.idx = c()
    tdm.data =c()
    tdm.data <- as.data.frame(tdm[i])
    train.idx <- sample(nrow(tdm.data), round(nrow(tdm.data)*0.8))
    test.idx<-(1:nrow(tdm.data))[-train.idx]
    
    if (i ==1){
      train.kknn <-tdm.data[train.idx,]
      test.kknn <- tdm.data[test.idx,]
    } else {
      train.kknn <- rbind(train.kknn,tdm.data[train.idx,])
      test.kknn <- rbind(test.kknn,tdm.data[test.idx,])}
  }
  
  train.kknn.nl <- train.kknn[,!colnames(train.kknn)%in%"targetCandidate"]
  
  test.cand<-test.kknn$targetCandidate
  test.cand.nl<-test.kknn[,!colnames(test.kknn)%in%"targetCandidate"]

  
  kknn.5<-kknn.ordinal(targetCandidate~., train.kknn,test.cand.nl,distance =2,k=7,kernel = "triangular",param=0.5)
  fit.test5 <- fitted(kknn.5)
  conf.mat5<-table(Actual=test.cand,"Predictions"= fit.test5)
  fit_one <- model.matrix(~0+fit.test5)
  colnames(fit_one) <- substr(colnames(fit_one),10,10)
  
  mauc_s$ori<-mauc(test.cand,fit_one)$mauc
  Gmetric_s$ori<-Gmetric(conf.mat5)
  dan_s$ori<-danger(train.kknn,train.kknn.nl)
  
  ##SMOTE FOR ORIGIN.DATA
  smote.kknn = c()
  maxnum <- summary(train.kknn$targetCandidate)
  maxnum <- as.numeric(maxnum[which.max(maxnum)])
  datasmote <- train.kknn
  datasmote$targetCandidate <- factor(datasmote$targetCandidate,ordered = FALSE )
  c(lv1,lv2,lv3,lv4,lv5) %<-% c(as.numeric(summary(train.kknn$targetCandidate)))
  C.perc = list("1" = maxnum/lv1,"2" = maxnum/lv2,"3" = maxnum/lv3,"4" = maxnum/lv4,"5" = maxnum/lv5)
  smote.kknn <- SmoteClassif(targetCandidate~., datasmote, C.perc, k=5)
  smote.kknn$targetCandidate <- as.ordered(smote.kknn$targetCandidate)
  smote.kknn.nl <- smote.kknn[,!colnames(smote.kknn)%in%"targetCandidate"]
  smote <- kknn.ordinal(targetCandidate~., smote.kknn,test.cand.nl,distance =2,k=5,kernel = "triangular",param=0.5)
  smote5 <- fitted(smote)
  smote.5<- table(Actual=test.cand,"Predictions"= smote5)
  fit_one <- model.matrix(~0+smote5)
  colnames(fit_one) <- substr(colnames(fit_one),10,10)
  
  mauc_s$smo<-mauc(test.cand,fit_one)$mauc
  Gmetric_s$smo<-Gmetric(smote.5)
  dan_s$smo<-danger(smote.kknn,smote.kknn.nl)
  
  
  #SCUT 
  scout.kknn = c()
  m <-  mean(summary(train.kknn$targetCandidate)) #2,3  <- under 1,4,5 <- over
  
  clus2.kknn <- createDataPartition(Mclust(subset(train.kknn,targetCandidate==2)[,-4])$classification,p=m/lv2)
  clus3.kknn <- createDataPartition(Mclust(subset(train.kknn,targetCandidate==3)[,-4])$classification,p=m/lv3)
  
  scout2.kknn <- subset(train.kknn,targetCandidate==2)[c(clus2.kknn$Resample1),]
  scout3.kknn <- subset(train.kknn,targetCandidate==3)[c(clus3.kknn$Resample1),]
  
  smote_s <- train.kknn
  smote_s$targetCandidate <- factor(smote_s$targetCandidate,ordered=FALSE)
  sm1.kknn <- SmoteClassif(targetCandidate~., subset(smote_s,targetCandidate==1), list("1" = m/lv1), k=7)
  sm4.kknn <- SmoteClassif(targetCandidate~., subset(smote_s,targetCandidate==4), list("4" = m/lv4), k=7)
  sm5.kknn <- SmoteClassif(targetCandidate~., subset(smote_s,targetCandidate==5), list("5" = m/lv5), k=7)
  
  scout.kknn <- rbind(scout2.kknn,scout3.kknn,sm1.kknn,sm4.kknn,sm5.kknn)
  scout.kknn$targetCandidate <- as.ordered(scout.kknn$targetCandidate)
  scout.kknn.nl <- scout.kknn[,!colnames(scout.kknn)%in%"targetCandidate"]
  
  scout <- kknn.ordinal(targetCandidate~., scout.kknn,test.cand.nl,distance =2,k=7,kernel = "triangular",param=0.5) 
  scout5 <- fitted(scout)
  scout.5<- table(Actual=test.cand,"Predictions"= scout5)
  fit_one <- model.matrix(~0+scout5)
  colnames(fit_one) <- substr(colnames(fit_one),10,10)
  
  mauc_s$scut<-mauc(test.cand,fit_one)$mauc
  Gmetric_s$scut<-Gmetric(scout.5)
  dan_s$scut<-danger(scout.kknn,scout.kknn.nl)

  
  #Wilson's editing . multiclass version
  ENNtrain <- train.kknn
  ENNtrain$targetCandidate <- factor(ENNtrain$targetCandidate,ordered=FALSE)
  #str(ENNtrain)
  ENN <- ENNClassif(targetCandidate~.,ENNtrain)
  edit.kknn <- ENN[[1]]
  edit.kknn$targetCandidate <- as.ordered(edit.kknn$targetCandidate)
  
  edit <- kknn.ordinal(targetCandidate~., edit.kknn,test.cand.nl,distance =2,k=7,kernel = "triangular",param=0.5)
  edit5 <- fitted(edit)
  edit.5<- table(Actual=test.cand,"Predictions"= edit5)
  fit_one <- model.matrix(~0+edit5)
  colnames(fit_one) <- substr(colnames(fit_one),10,10)
  
  mauc_s$edit<-mauc(test.cand,fit_one)$mauc
  Gmetric_s$edit<-Gmetric(edit.5)
  dan_s$edit<-danger(edit.kknn,edit.kknn[,-4])
  
  
  #InCuBAtE
  

  weights <- relief(targetCandidate~., data.frame(train.kknn), neighbours.count = 5, sample.size = 10)
  features <- cutoff.k.percent(weights, 0.7) 
  
  
  In_test <- test.cand.nl[,c(features)]
  In_train <- train.kknn.nl[,c(features)]
  In_train <- cbind(In_train, train.kknn$targetCandidate)
  colnames(In_train)[dim(In_train)[2]] <- "targetCandidate"
  
  maxnum <- summary(train.kknn$targetCandidate)
  maxnum <- as.numeric(maxnum[which.max(maxnum)])
  
  du_train <- rbind(In_train,In_train)
  class <- as.numeric(unique(train.kknn$targetCandidate))
  class <- class[-which(summary(train.kknn$targetCandidate)==maxnum)]
  
  for(kk in class){
    n <- dim(subset(du_train,du_train$targetCandidate==kk))[1]
    N <- dim(subset(train.kknn,train.kknn$targetCandidate==kk))[1]
    
    while((n-N) <maxnum){
      maximum <- as.numeric(apply(subset(du_train,du_train$targetCandidate==kk),2,max))
      minimum <- as.numeric(apply(subset(du_train,du_train$targetCandidate==kk),2,min))
      g_data <- c() 
      for(s in 1:(dim(In_train)[2]-1)){
        if(minimum[s] > maximum[s]){
          minimum[s] <- min(subset(du_train[,s],du_train$targetCandidate==kk))
          maximum[s] <- max(subset(du_train[,s],du_train$targetCandidate==kk))
        }
        if(minimum[s]==maximum[s]){
          g_data[s] <- minimum[s]
        } else {g_data[s] <- runif(1,min=minimum[s],max=maximum[s])} 
      }
      g_data[dim(In_train)[2]] <- kk
      
      syn_data <- rbind(du_train,g_data)
      du_train$targetCandidate <- as.ordered(du_train$targetCandidate)
      prediction <- kknn.ordinal(targetCandidate~.,du_train,syn_data[dim(syn_data)[1],-dim(syn_data)[2]],distance =2,k=7,kernel = "triangular",param=0.5)$fitted.values 
      
      if(prediction==kk){
        du_train <- syn_data
      }
      
      n <- dim(subset(du_train,du_train$targetCandidate==kk))[1]
    }
    
  }
  
  Incu_train <- du_train[-c(1:nrow(In_train)),]
  
  IncuBAtE <- kknn.ordinal(targetCandidate~., Incu_train,In_test,distance =2,k=7,kernel = "triangular",param=0.5)
  IncuBAtE5 <- fitted(IncuBAtE)
  IncuBAtE.5<- table(Actual=test.cand,"Predictions"= IncuBAtE5)
  fit_one <- model.matrix(~0+IncuBAtE5)
  colnames(fit_one) <- substr(colnames(fit_one),10,10)
  
  mauc_s$IncuBAtE<-mauc(test.cand,fit_one)$mauc
  Gmetric_s$IncuBAtE<-Gmetric(IncuBAtE.5)
  dan_s$IncuBAtE<-danger(Incu_train,Incu_train[,-dim(Incu_train)[2]])
  
  #SVD
  classattribute <- 4
  classnumber <- 5
  
  
  train_kknn<-train.kknn
  train_kknn$targetCandidate<-as.ordered(as.numeric(train_kknn$targetCandidate))
  
  
  listdd<-list()
  groups<-createFolds(train_kknn[,classattribute], k = 5, list = TRUE, returnTrain = FALSE)
  for(rr in 1:5){
    betatrain <- train_kknn[-groups[[rr]],]
    betavalidation <- train_kknn[groups[[rr]],]
    listdd[[rr]]<-list(betatrain=betatrain,betavalidation=betavalidation)
  }
 
  fitness<-function(par,listdd){
    ss<-0
    for(aa in 1:5){
      ss<-ss+soft_thresholding(par,listdd[[aa]])
    }
    return(ss)
  }
  
  GA1 <- ga(type = "real-valued", 
            fitness =  fitness, listdd=listdd, optim = TRUE, 
            lower = c(9,9,9,9,9), upper = c(50,50,50,50,50), 
            popSize = 50,maxiter = 100, run = 30)    
  
  
  
  bound <- summary(GA1)$solution[1,]

  beta_bound <- rbind(beta_bound,data.frame(matrix(bound,ncol=classnumber,byrow = T)))
  
  
  
  syntrain.kknn <- c()
  syntrain.kknn <- train.kknn
  syntrain.kknn[,-4] <-SVDtruncation(bound)
  syntrain.kknn.nl <- syntrain.kknn[,!colnames(syntrain.kknn)%in%"targetCandidate"]

  plot(syntrain.kknn[,2],syntrain.kknn[,3],col=(as.numeric(syntrain.kknn$targetCandidate)))
  plot(train.kknn[,2],train.kknn[,3],col=(as.numeric(train.kknn$targetCandidate)))
  

  
  
  syn.kknn <-kknn.ordinal(targetCandidate~., syntrain.kknn,test.cand.nl,distance =2,k=7,kernel = "triangular",param=0.5)
  fit.syn <- fitted(syn.kknn)
  conf.syn<-table(Actual=test.cand,"Predictions"= fit.syn)
  fit_one <- model.matrix(~0+ fit.syn )
  colnames(fit_one) <- substr(colnames(fit_one),10,10)
  
  mauc_s$syn<-mauc(test.cand,fit_one)$mauc
  Gmetric_s$syn<-Gmetric(conf.syn)
  dan_s$syn<-danger(syntrain.kknn,syntrain.kknn.nl)
  

  bound3 <- tuningbeta_whole(train.kknn,4)
  beta_bound3 <- c(beta_bound3,bound3)
  
  syntrain3 <- c()
  syntrain3 <- rbind(train.kknn,test.kknn)
  syntrain3[,-4] <-dataapproximation(rbind(train.kknn,test.kknn),bound3)
  syntrain.kknn3 <- syntrain3[1:nrow(train.kknn),]
  syntest.kknn3 <- syntrain3[(nrow(train.kknn)+1):nrow(syntrain3),]
  syntrain.kknn.nl3 <- syntrain.kknn3[,!colnames(syntrain.kknn3)%in%"targetCandidate"]
  
  syn.kknn <-kknn.ordinal(targetCandidate~., syntrain.kknn3,test.cand.nl,distance =2,k=7,kernel = "triangular",param=0.5)
  fit.syn <- fitted(syn.kknn)
  conf.syn<-table(Actual=test.cand,"Predictions"= fit.syn)
  fit_one <- model.matrix(~0+ fit.syn )
  colnames(fit_one) <- substr(colnames(fit_one),10,10)
  if ((dim(fit_one)[2] < 3)){
    mauc_s$syn3<-0
    Gmetric_s$syn3<-0
  } else {
    mauc_s$syn3<-mauc(test.cand,fit_one)$mauc
    Gmetric_s$syn3<-Gmetric(conf.syn)
  }
  
  dan_s$syn3<-danger(syntrain.kknn3,syntrain.kknn.nl3)

  maucs<-rbind(maucs,data.frame(mauc_s))
  Gmetrics<-rbind(Gmetrics,data.frame(Gmetric_s))
  dans<-rbind(dans,data.frame(dan_s))
}


obj<-list(maucs,Gmetrics,dans,beta_bound,beta_bound3)
round(colMeans(maucs),3)
round(colMeans(Gmetrics),3)
round(colMeans(dans),3)
round(colMeans(beta_bound),3)

round(apply(maucs,2,sd),3)
round(apply(Gmetrics,2,sd),3)
round(apply(dans,2,sd),3)


save(obj,file="real_1022.RData")

mauc<-obj[[1]]
gmean<-obj[[2]]
dan<-obj[[3]]

boxplot(mauc[,c(1:7)],names=c("Original","SVSMC","SVD","SMOTE","SCUT","ENN","InCuBAtE"))
boxplot(gmean[,c(1:7)],names=c("Original","SVSMC","SVD","SMOTE","SCUT","ENN","InCuBAtE"))
boxplot(dan[,c(1:7)],names=c("Original","SVSMC","SVD","SMOTE","SCUT","ENN","InCuBAtE"))

boxplot(beta_bound)


## visualization
library(wordcloud)
for(ii in 1:5){
  level_1<-tdm.stack[tdm.stack[,4]==ii,-4]
  d <- data.frame(word = names(level_1),freq=as.numeric(apply(level_1,2,sum)))
  #png(paste("map",ii,".png",sep=""), width = 600, height = 600)
  wordcloud(words = d$word, freq = d$freq,scale=c(8,.5),min.freq = 0,
            max.words=100, random.order=FALSE, rot.per=0.35, 
            colors=brewer.pal(8, "Dark2"))
  #dev.off()
}

hh<-data.frame(var=c(rep(1,161),rep(2,397),rep(3,209),rep(4,94),rep(5,27)))
library(ggplot2)
ggplot(hh,aes(x= var))+geom_histogram(bins=5,color="white")+xlab("Level")+ylab("Number of events")+theme(text = element_text(size=15))



##

DT<-(data.frame(train.kknn[,c(2,17,4)]))
DTu<-ddply(DT,.(area,killed,targetCandidate),nrow)
names(DTu)[3:4]<-c("Level","Weight")
DTu$Level<-as.character(DTu$Level)

p <- ggplot(DTu, aes(y=area, x=killed))
p+ geom_point(aes(size = Weight,shape=Level,color=Level), alpha = 0.8) +
  scale_fill_viridis_d(option = "inferno") +
  scale_size_continuous(range = c(1, 20))

##
DT2<-(data.frame(syntrain.kknn[,c(2,17,4)]))
DTu2<-ddply(DT2,.(area,killed,targetCandidate),nrow)
names(DTu2)[3:4]<-c("Level","Weight")
DTu2$Level<-as.character(DTu2$Level)

p2 <- ggplot(DTu2, aes(y=area, x=killed))
p2+ geom_point(aes(size = Weight,shape=Level,color=Level), alpha = 0.7) +
  scale_fill_viridis_d(option = "inferno") +
  scale_size_continuous(range = c(1, 20))

