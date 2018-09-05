# Machine Learning scripts
# by: Ranko Gacesa, UMCG (2017)
# ==============================================

library(caret)
library(plotROC)
library(pROC)

# do train & testing, make confusion matrices

doTrainTest <- function(dataI,mdl,trainPerc=0.7) {
  inTrain <- createDataPartition(y=dataI$Diagnosis,p=trainPerc,list=F)
  trainSet <-dataI[inTrain,]
  testSet <- dataI[-inTrain,]
  t <- train(Diagnosis ~ ., data=trainSet, method=mdl)
  cTest <- confusionMatrix(predict(t,testSet),testSet$Diagnosis)
  cTrain <- confusionMatrix(predict(t,trainSet),trainSet$Diagnosis)
  return (list(t,cTrain,cTest))
}


# ============================
#
# plot all learning curves
#
# ===========================
plotAllLearningCurves <- function(dataI,dataName="",dataNameShort="",
                                  mdls=c("glm","gbm","rpart2","avNNet","nnet","pcaNNet","hdrda","svmRadial","svmLinear3",
                                         "rda","rrlda","regLogistic","rf","RRFglobal"),
                                  bts=3,sStep=2,sMin=50,outTable=T,outPlots=T,outFolder='plots',trainRepNR = 2,trainBootNR = 10,
                                  doROC = T,posClass="IBS") {
  # outtable
  outTbl = data.frame(Model=character(),Metric=numeric(),Value=numeric(),SD=numeric())
  # go over all models
  for (md in mdls) {
    startTime <- Sys.time()
    print (paste(">>>>>> PREPPING LEARNING CURVE FOR ",md," <<<<<<"))
    #registerDoSNOW(makeCluster(4, type = "SOCK"))
    lCurve <- prepLearningCurve(dataIn = dataI,mdl = md, minSam = sMin,boots = bts,samStep = sStep,trNumber = trainRepNR,trBoot = trainBootNR,
                                saveVarImp=paste(outFolder,'/lCurve_',dataNameShort,'_',md,'_varImp.png',sep=''),
                                saveVarImpTit=paste('Covariate importance (',md,')',sep=''),
                                posClass=posClass,ROCtitle=paste('ROC: ',dataName,' [',md,']',sep='') )
    g <- plotLearningCurves(lCurve,tit=paste("L-Curve",dataName,md))
    ggsave(g,filename = paste(outFolder,'/lCurve_',dataNameShort,'_',md,'.png',sep=''),width = 9,height = 9)
    endTime <- Sys.time()
    fStats <- as.data.frame(lCurve[nrow(lCurve),3:ncol(lCurve)])
    fStats$TIME <- round(as.numeric(difftime(time1 = endTime, time2= startTime,units = "sec")))
    fStats$TIME.SD <- 0.0
    
    fStats$Model <- md
    fStatsSD <- fStats[,grep("\\.SD",colnames(fStats))]
    fStatsSD$Model <- md
    fStats <- fStats[,grep("\\.SD",colnames(fStats),invert = T)]
    fStatsLong <- gather(fStats,key = "Metric",value = "Value",Acc:TIME)
    fStatsSDLong <- gather(fStatsSD,key = "Metric",value = "SD",Acc.SD:TIME.SD)
    t <- cbind.data.frame(fStatsLong,fStatsSDLong$SD)
    t[is.na(t)] <- 0.0
    colnames(t) <- c("Model","Metric","Value","SD")
    
    outTbl <- rbind.data.frame(outTbl,t)
    
    #outTbl$Model <- as.character(outTbl$Model)
    #print(outTbl)
  }
  #for (i in c(2:9)) {outTbl[[i]] <- as.numeric(outTbl[[i]])}
  
  g <- ggplot(data=outTbl[outTbl$Metric!="TIME",],aes(x=reorder(Model,Value),y=Value,col=reorder(Model,Value))) + geom_point(size=3) + facet_grid(Metric ~ .) + 
    geom_errorbar(aes(ymin=Value-SD, ymax=Value+SD), width=0.15) + ylim(0,1) + ggtitle(paste(dataName," Metrics",sep="")) + 
    xlab("ML algorithm") + ylab("Prediction metric value") + guides(col=guide_legend(title="ML algorithm")) + theme(axis.text.x = element_text(size = rel(0.75), angle = 00))
  
  ggsave(g,filename = paste(outFolder,'/lCurve_',dataNameShort,'__COMP.png',sep=''),width = 12,height = 10)
  
  g <- ggplot(data=outTbl[outTbl$Metric=="TIME",],aes(x=reorder(Model,Value),y=Value,fill=reorder(Model,Value))) + geom_col() + 
    ggtitle(paste(dataName," Training Time",sep="")) + xlab("ML algorithm") + ylab("Training time (seconds)") + 
    guides(fill=guide_legend(title="ML algorithm")) + theme(axis.text.x = element_text(size = rel(0.75), angle = 00))
  ggsave(g,filename = paste(outFolder,'/lCurve_',dataNameShort,'__TIME.png',sep=''),width = 12,height = 8)
  
  
  if (outTable) {
    write.csv(outTbl,file=paste(outFolder,'/lCurve_',dataNameShort,'__COMP_TABLE.csv',sep=''))
  }

  if (outTable) {
    return(outTbl)
  }
}

# ===============================================================================
# generate ROC based on RF trained on given model (70% of data)
#
# ===============================================================================
generateROC <- function(dataModel,trSet=0.7,klasses=c("IBS","HC"),tit="") {
  inTrain <- createDataPartition(y=dataModel$Diagnosis,p=trSet,list=F)
  trainSet <-dataModel[inTrain,]
  testSet <- dataModel[-inTrain,]
  fittedModel <- train(Diagnosis ~ ., data=trainSet,method="rf",metric="ROC",
                trControl=trainControl(savePrediction=T,classProbs=T,summaryFunction=twoClassSummary))
  sI <- fittedModel$pred$mtry == fittedModel$bestTune$mtry
  
  grug <- fittedModel$pred[sI,grep(klasses[1],colnames(fittedModel$pred))]
  g <- ggplot(fittedModel$pred[sI,], aes(m=grug, d=factor(obs, levels = klasses))) + 
    style_roc(xlab = "1 - Specificity / False Positive rate",
              ylab = "Sensitivity / True Positive rate") + 
    geom_roc(n.cuts = 0,color="blue",size=1.5) + 
    coord_equal() 
  auc <- fittedModel$results[fittedModel$results$mtry==fittedModel$bestTune$mtry,]$ROC
  g <- g + annotate("text", x=0.70, y=0.3, label=paste("AUC =", round(auc, 3)))
  g <- g + ggtitle(tit)
  g
}

# ===============================================================================
# generate ROC for listed model (mdl) and testSet (testSet)
# positive class is inputed one (posClass), other class is considered to be negative
# ===============================================================================
generateMdlROC <- function(mdl,testSet,posClass="IBS",tit="") {
  predP <- predict(mdl,testSet,type="prob")
  r <- roc(testSet$Diagnosis,predP[[posClass]],auc=T,percent=F)
  #cir <- ci(r,of="se", sp = seq(0, 100, 5), boot.n=100)
  #plot(r,print.auc=T,col="darkgreen")
  
  g <- pROC::ggroc(r,col="darkblue",size=1.5,alpha=0.75) +
  annotate("text", x = .75, y = .25, label = paste("AUC =",round(r$auc[1],2) )) +
    ggtitle(tit)

}
# ===============================================================================
# generate ROC for listed model (mdl) and testSet (testSet) + trainSet (trainSet)
# positive class is inputed one (posClass), other class is considered to be negative
# ===============================================================================
generateMdlROCTrainTest <- function(mdl,trainSet,testSet,posClass="IBS",tit="",doSmooth=F) {
  
  predTr <- predict(mdl,trainSet,type="prob")
  rTr <- roc(trainSet$Diagnosis,predTr[[posClass]],auc=T,percent=F)
  
  predTst <- predict(mdl,testSet,type="prob")
  rTst <- roc(testSet$Diagnosis,predTst[[posClass]],auc=T,percent=F)
  
  #g <- pROC::ggroc(r,col="darkblue",size=1.5,alpha=0.75) +
  #  annotate("text", x = .75, y = .25, label = paste("AUC =",round(r$auc[1],2) )) +
  #  ggtitle(tit)  
  if (doSmooth) {
    g <- pROC::ggroc(list(Train=rTr,Test=rTst), size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(list(Train=rTr,Test=rTst), size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(tit) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    annotate("text", x = .45, y = .23, label = paste("Training Set AUC =",round(rTr$auc[1],2)),hjust = 0) +
    annotate("text", x = .45, y = .30, label = paste("Test Set AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    geom_abline(slope = 1,intercept = 1,col="grey") +
    ggtitle(tit) + xlab("Specificity") + ylab("Sensitivity") +
    labs(color='Dataset')  
  g
}
# ===============================================================================
# generator of comparative ROC for MULTIPLE MODELS and ONE DATASET (test set)
# ===============================================================================
generateMdlROCComparative <- function(mdls,modelNames,testSet,posClass="IBS",tit="",doSmooth=F) {
  rocs = list()
  c = 0
  for (m in mdls) {
    c = c + 1
    namen = modelNames[c]
    pr <- predict(m,testSet,type="prob")
    rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
  }
  
  if (doSmooth) {
  g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(tit) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    geom_abline(slope = 1,intercept = 1,col="grey") +
    labs(color='Dataset') 
  # annotate w auc
  c = 0.0
  for (i in order(modelNames,decreasing = T)) {
    c = c + 1.0
    g <- g + annotate("text", x = .35, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0)
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  g
}
# ===============================================================================
# generator of comparative ROC for MULTIPLE DATASETS and ONE MODEL
# ===============================================================================
generateMdlROCCompareData <- function(mdl,testSetNames,testSets,posClass="IBS",tit="",doSmooth=F) {
  rocs = list()
  c = 0
  m <- mdl
  for (testSet in testSets) {
    c = c + 1
    namen = testSetNames[c]
    pr <- predict(m,testSet,type="prob")
    rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
  }
  
  if (doSmooth) {
    g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(tit) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    geom_abline(slope = 1,intercept = 1,col="grey") +
    labs(color='Dataset') 
  # annotate w auc
  c = 0.0
  for (i in order(testSetNames,decreasing = T)) {
    c = c + 1.0
    g <- g + annotate("text", x = .35, y = 0.00+0.06*c, label = paste(testSetNames[i],"AUC =",round(rocs[[testSetNames[i] ]]$auc[1],2)),hjust = 0)
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  g
}

# ===============================================================================
# compares model performance on different datasets
# ===============================================================================
compareMdlDatasets <- function(mdl,dataSets,dataSetNames,posClass="IBS",roc.smooth=F,tit="",
                               response="Diagnosis",roc.conf = T,roc.conf.boot = 100,target = F) {
  trainedModels <- list()
  c = 0
  # list of ROC curves
  rocs = list()
  # data frame of predictions
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  # generate results for each dataset
  for (ds in dataSets) {
    c = c + 1
    testSet <- ds
    fittedModel <- mdl
    trainedModels[[c]] <- fittedModel
    pr <- predict(fittedModel,newdata = testSet,type="prob")
    namen = dataSetNames[c]
    if (roc.smooth) {
      rocs[[namen]] <- smooth(roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F))
    } else {
      rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
    }
    pr2 <- predict(fittedModel,testSet)
    conf <- confusionMatrix(pr2,testSet[[response]])
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="F1","Value"=conf$byClass[["F1"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataSet"=namen,"Metric"="Prevalence","Value"=conf$byClass[["Prevalence"]],stringsAsFactors=F))
  }
  
  # extract coordinates
  if (!roc.conf) {
    rocsDFprecalc = NULL
    rocStep = 0.01
    for (i in seq(1,length(rocs))) {
      if (class(rocsDFprecalc) != "data.frame") {
        if (!roc.smooth) {
          rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalc) <- NULL
        } else {
          rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x=c(seq(0,0.1,rocStep/10),seq(0.1,1,rocStep)),ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalc) <- NULL
        }
        rocsDFprecalc$Dataset <- names(rocs)[i]
        rocsDFprecalc$sz = 1
      } else {
        if (!roc.smooth) {
          rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalctmp) <- NULL
        } else {
          rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x=c(seq(0,0.1,rocStep/10),seq(0.1,1,rocStep)),ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalctmp) <- NULL
        }
        rocsDFprecalctmp$Dataset <- names(rocs)[i]
        rocsDFprecalctmp$sz = 1
        rocsDFprecalc <- rbind(rocsDFprecalc, rocsDFprecalctmp)
        rownames(rocsDFprecalc) <- NULL
      }
    }
  }
  
  if (roc.conf) {
    rocDFConf <- F
    rocStep = 0.01
    boots = roc.conf.boot
    meanConf <- F
    for (i in seq(1,length(rocs))) {
      if (roc.smooth) {
        sens.ci <- ci.sp(smooth(rocs[[i]]), sensitivities =c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep)),
                         conf.level=0.66,boot.n = boots)
      } else {
        sens.ci <- ci.sp((rocs[[i]]), sensitivities=c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep)),
                         conf.level=0.66,boot.n = boots)
      }
      rocConf <- as.data.frame(sens.ci)
      rocConf$sensitivity <- as.numeric(rownames(rocConf))
      colnames(rocConf) <- c("sp.low","specificity","sp.high","sensitivity")
      rocConf$Dataset <- names(rocs)[i]
      rocConf$sz = 1.0
      rownames(rocConf) <- NULL
      if (class(meanConf) != "data.frame") {
        meanConf <- rocConf[,c(4,1,2,3)]
      } else {
        meanConf <- cbind(meanConf,rocConf[,c(1,2,3)])
      }
      if (class(rocDFConf) != "data.frame") {
        rocDFConf <- rocConf
      } else {
        rocDFConf <- rbind.data.frame(rocDFConf,rocConf)
      }
      if (!roc.smooth) {
        rocDFConf <- rbind.data.frame(rocDFConf,data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,Dataset=names(rocs)[i],sz=1))
      }
    }
    rocsDFprecalc <- rocDFConf
  }
  
  if (roc.conf) {
    rocsDFprecalc <- rocsDFprecalc[order(rocsDFprecalc$specificity,decreasing = F),]
    g <- ggplot(rocsDFprecalc,aes(y=specificity,x=sensitivity,col=Dataset)) +
      geom_ribbon(data=rocsDFprecalc,aes(ymin=sp.low,ymax=sp.high,fill=Dataset),alpha=0.2,colour=NA) + geom_line(size=1.25)
  } else {
    rocsDFprecalc <- rocsDFprecalc[order(rocsDFprecalc$sensitivity,decreasing = F),]
    g <- ggplot(rocsDFprecalc,aes(y=specificity,x=sensitivity,col=Dataset)) + geom_line(size=1.25) #pROC::ggroc(rocs) + geom_line(size=1.25)
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1) +
    ggtitle(paste(tit,"ROC") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +  
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color='Dataset') 
  # annotate w auc
  c = 0.0
  for (i in order(dataSetNames,decreasing = T)) {
    c = c + 1.0
    g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(dataSetNames[i],"AUC =",round(rocs[[dataSetNames[i] ]]$auc[1],2)),hjust = 0)
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  # targets
  tarDF <- F  
  if (target) {
    for (u in unique(predDF$DataSet )) {
      if (u != "multi.class") {
        if (class(tarDF) != "data.frame") {
          tarDF <- data.frame(DataSet=u,
                              Sensitivity=predDF[predDF$DataSet==u & predDF$Metric=="Sensitivity",]$Value,
                              Specificity=predDF[predDF$DataSet==u & predDF$Metric=="Specificity",]$Value,
                              F1=predDF[predDF$DataSet==u & predDF$Metric=="F1",]$Value,
                              B.Acc=predDF[predDF$DataSet==u & predDF$Metric=="B.Acc",]$Value
          )
        } else {
          tarDF <- rbind.data.frame(tarDF,data.frame(DataSet=u,
                                                     Sensitivity=predDF[predDF$DataSet==u & predDF$Metric=="Sensitivity",]$Value,
                                                     Specificity=predDF[predDF$DataSet==u & predDF$Metric=="Specificity",]$Value,
                                                     F1=predDF[predDF$DataSet==u & predDF$Metric=="F1",]$Value,
                                                     B.Acc=predDF[predDF$DataSet==u & predDF$Metric=="B.Acc",]$Value))
        }
      }
    }
    g <- g + geom_point(data=tarDF,aes(x=Specificity,y=Sensitivity,col=DataSet),size=4,shape=4,stroke=2)
  }
  
  # do metrics plot
  unik <- length(dataSets)
  gg2 <- ggplot(data=predDF,aes(col=Metric,y=Value,x=DataSet,shape=Metric)) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + 
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"Prediction metrics") ) +
    geom_vline(xintercept =  c(1:length(dataSetNames)) ,size=max(35,130-unik*20),alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  list(g,gg2,predDF)
}


# ===============================================================================
# compares training CV of multiple trained models
# -> dataMdls should be list( dataframes(Diagnosis ~ covariate) )
# -> note: models should predict Diagnosis, with positive class = posClass
# ===============================================================================
compareModelsTrainingCV <- function(fittedMdls,modelNames,mtd="glm",posClass="IBS",doSmooth=F,tit="",annotateAUConly=F) {
  c = 0
  rocs = list()
  for (m in fittedMdls) {
    c = c + 1
    namen = modelNames[c]
    rocs[[namen]] <- roc(predictor = fittedMdls[[c]]$pred[[posClass]], response = fittedMdls[[c]]$pred$obs,auc=T,percent = F)
  }
  
  # do combined ROC plot
  if (doSmooth) {
    g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(paste(tit,"ROC (",mtd,")") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color='Dataset') 
  # annotate w auc
  c = -1
  for (i in order(modelNames,decreasing = T)) {
    c = c + 1.0
    if (!annotateAUConly) {
      g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(modelNames[i]," AUC = ",round(rocs[[modelNames[i] ]]$auc[1],2),
                                                                       "; Kappa = ",round(getTrainPerf(fittedMdls[[i]])[[2]],2),
                                                                       "; ACC = ",round(getTrainPerf(fittedMdls[[i]])[[1]],2),sep=""),hjust = 0)
    } else {
      g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0)
      
    }
#    c = c + 1.0
#    g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0)
  }
  #return plot
  g
}

# ===============================================================================
# compares multiple data models - each one is built from its own dataset
# does the model training
# -> dataMdls should be list( dataframes(Diagnosis ~ covariate) )
# -> note: models should predict Diagnosis, with positive class = posClass
# ===============================================================================
trainCompareModels <- function(dataMdls,modelNames,mtd="glm",posClass="IBS",doSmooth=F,trSet=0.7,tit="") {
  c = 0
  rocs = list()
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  #trainedModels = list()
  for (m in dataMdls) {
    c = c + 1
    inTrain <- createDataPartition(y=m$Diagnosis,p=trSet,list=F)
    trainSet <- m[inTrain,]
    testSet <- m[-inTrain,]
    fittedModel <- train(Diagnosis ~ ., data=trainSet,method=mtd,metric="Kappa")
    pr <- predict(fittedModel,testSet,type="prob")
    namen = modelNames[c]
    rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
    pr2 <- predict(fittedModel,testSet)
    conf <- confusionMatrix(pr2,testSet$Diagnosis)
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Det.Rate","Value"=conf$byClass[["Detection Rate"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Prevalence","Value"=conf$byClass[["Prevalence"]],stringsAsFactors=F))
  }
  # do combined ROC plot
  if (doSmooth) {
    g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(paste(tit,"ROC (",mtd,")") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color='Data Model') 
  # annotate w auc
  c = 0.0
  for (i in order(modelNames,decreasing = T)) {
    c = c + 1.0
    g <- g + annotate("text", x = .35, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0)
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  # do metrics plot
  gg2 <- ggplot(data=predDF,aes(col=Metric,y=Value,x=DataModel,shape=Metric)) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + xlab("Data Model") +
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"Prediction metrics (",mtd,")") ) +
    geom_vline(xintercept =  c(1:length(modelNames)) ,size=40,alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  list(g,gg2,predDF)
}




# ===============================================================================
# compares multiple models on multiple test datasets
# (each model with paired with appropriate test-set in order they are entered)
# does not do model training
# -> fittedMdls should be list( fitted models ) generated by caret train
# ===============================================================================
compareModelsOnDatasets <- function(fittedMdls,testSets,modelNames,mtd="glm",posClass="IBS",doSmooth=F,tit="",annotateAUConly=F) {
  c = 0
  rocs = list()
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  for (m in fittedMdls) {
    c = c + 1
    #trainedModels[[c]] <- m
    pr <- predict(m,testSets[[c]],type="prob")
    namen = modelNames[c]
    rocs[[namen]] <- roc(testSets[[c]]$Diagnosis,pr[[posClass]],auc=T,percent=F)
    pr2 <- predict(m,testSets[[c]])
    conf <- confusionMatrix(pr2,testSets[[c]]$Diagnosis)
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Det.Rate","Value"=conf$byClass[["Detection Rate"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Prevalence","Value"=conf$byClass[["Prevalence"]],stringsAsFactors=F))
  }
  
  # do combined ROC plot
  if (doSmooth) {
    g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(paste(tit,"(",mtd,")") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color='Dataset') 
  # annotate w auc
  c = -1.0
  for (i in order(modelNames,decreasing = T)) {
    c = c + 1.0
    if (annotateAUConly) {
      g <- g + annotate("text", x = .35, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0)
    } else {
      g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(modelNames[i]," AUC = ",round(rocs[[modelNames[i] ]]$auc[1],2),
                                                                     "; Kappa = ",round(predDF[predDF$DataModel==modelNames[i] & predDF$Metric=="Kappa",]$Value,2),
                                                                     "; ACC = ",round(predDF[predDF$DataModel==modelNames[i] & predDF$Metric=="Acc",]$Value,1),sep=""),hjust = 0)
    }
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  # do metrics plot
  gg2 <- ggplot(data=predDF,aes(col=Metric,y=Value,x=DataModel,shape=Metric)) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + 
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"(",mtd,")") ) +
    geom_vline(xintercept =  c(1:length(modelNames)) ,size=40,alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  list(g,gg2,predDF)
}



# ===============================================================================
# compares multiple models on one dataset
# does not do model training
# -> fittedMdls should be list( fitted models ) generated by caret train
# ===============================================================================
compareModelsOnDataset <- function(fittedMdls,testSet,modelNames,mtd="glm",posClass="IBS",doSmooth=F,tit="") {
  #trainedModels <- list()
  c = 0
  rocs = list()
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  for (m in fittedMdls) {
    c = c + 1
    #trainedModels[[c]] <- m
    pr <- predict(m,testSet,type="prob")
    namen = modelNames[c]
    rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
    pr2 <- predict(m,testSet)
    conf <- confusionMatrix(pr2,testSet$Diagnosis)
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Det.Rate","Value"=conf$byClass[["Detection Rate"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Prevalence","Value"=conf$byClass[["Prevalence"]],stringsAsFactors=F))
  }
  
  # do combined ROC plot
  if (doSmooth) {
    g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(paste(tit,"ROC (",mtd,")") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color='Dataset') 
  # annotate w auc
  c = 0.0
  for (i in order(modelNames,decreasing = T)) {
    c = c + 1.0
    g <- g + annotate("text", x = .35, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0)
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  # do metrics plot
  gg2 <- ggplot(data=predDF,aes(col=Metric,y=Value,x=DataModel,shape=Metric)) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + 
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"Prediction metrics (",mtd,")") ) +
    geom_vline(xintercept =  c(1:length(modelNames)) ,size=40,alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  list(g,gg2,predDF)
}

# ===============================================================================
# find correlated covariates
# note: input data should not have NAs
# ===============================================================================
findCorrelatedCovariates <- function(inData,cutOff) 
{
  numa <- sapply(inData, is.numeric)
  #print (numa)
  if (sum(numa) > 1) {
    correlationMatrix <- cor(inData[,sapply(inData, is.numeric)])
    cs <- colnames(inData)[sapply(inData, is.numeric)]
    # calculate correlation matrix
    #print(correlationMatrix)
    correlationMatrix[is.na(correlationMatrix)] <- 0.0
    # find attributes that are highly corrected (ideally >0.75)
    highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=cutOff,verbose = TRUE)
    # print indexes of highly correlated attributes
    print(cs[highlyCorrelated])
  } else {
    print ('less then 2 numeric covariates!')
  }
}

findCorrelatedCovariates2 <- function(inData,cutOff) 
{
  correlationMatrix <- cor(inData[,sapply(inData, is.numeric)])
  cs <- colnames(inData)[sapply(inData, is.numeric)]
  for (i in 1:length(cs)) {
    cs[i] <- strsplit(cs[i],'\\.')[[1]][length(strsplit(cs[i],'\\.')[[1]])]
    #print (cs[i])
  }
  # calculate correlation matrix
  #print(correlationMatrix)
  correlationMatrix[is.na(correlationMatrix)] <- 0.0
  cntr = 0
  for (i in 1:nrow(correlationMatrix)){
    for (j in 1:i) {
      if ( !is.na(correlationMatrix[i,j]) & !(i == j)){
        if ( correlationMatrix[i,j] > cutOff) {
          cntr = cntr + 1
          print(paste(cntr,cs[i], "-" , cs[j], ": ", correlationMatrix[i,j]))
        }
      }
    }
  }
}


#' Prepares learning curves (but does not draw them)
#'
#'
#'
prepLearningCurve <- function(dataIn,mdl,boots=1,trSet=0.7,samStep=1.5,minSam=25,scaleCenter=T,trBoot = 25,trNumber = 5,trainType="boot",
                              pProc = c("center","scale"), saveVarImp="", saveVarImpTit="",responseVar="Diagnosis",posClass="IBS",ROCtitle='')
{
  inTrain <- createDataPartition(y=dataIn[[responseVar]],p=trSet,list=F)
  trainSet <- dataIn[inTrain,]
  testSet <- dataIn[-inTrain,]
  
  p <- c()
  pp <- minSam / nrow(trainSet)
  if (pp <= 0) {pp = 0.01}
  while (pp <= 1.5) {
    p <- c(p,pp)
    pp <- pp * samStep
    if (pp > 1.0) {pp <- 1.0; p <- c(p,pp); break}
  }
  #print(p)
  #for (i in c(round(log(minSam,base=samStep)-0.5,digits=0):round(log(nrow(trainSet),base=samStep)+0.5,digits = 0)) ) {
  #  p <- c(p,min(samStep**i,nrow(trainSet)))
  #}
  
  results = data.frame()
  results <- rbind.data.frame(results,c(0.0,"Train",0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),stringsAsFactors = F)
#  results <- as.data.frame(apply(results,MARGIN = 2,FUN= function(x) as.character(x)))
  results <- rbind.data.frame(results,c(0.0,"Test",0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),stringsAsFactors = F)
#  results <- as.data.frame(apply(results,MARGIN = 2,FUN= function(x) as.character(x)))
  
  cnt = 0
  for (n in p) { # numbers loop
    #print(n)
    # vectors for accuracy measures
    cnt = cnt + 1
    acc = c()
    acctrain = c()
    kappa = c()
    kappatrain = c()
    sens <- c()
    spec <- c()
    f1 <- c()
    ppv <- c()
    npv <- c()
    senstrain <- c()
    spectrain <- c()
    f1train <- c()
    ppvtrain <- c()
    npvtrain <- c()
    auctrain <- c()
    auc <- c()
    
    for (b in c(1:boots)) { # bootstraps loop
      # another partition
      inTrain <- createDataPartition(y=dataIn[[responseVar]],p=trSet,list=F)
      trainSet <-dataIn[inTrain,]
      testSet <- dataIn[-inTrain,]
      # subselect
      #smplTr <- createDataPartition(y=trainSet[[responseVar]],p=trSet,list=F)
      #smpl <- trainSet[sample(nrow(trainSet), n), ]
      smplTake <- createDataPartition(y=trainSet[[responseVar]],p=n,list=F)
      smpl <- trainSet[smplTake,]
      
      print(paste("Test",cnt,"(",nrow(smpl),"cases);","bootstrap",b))
      # prep traincontrol
      tC = trainControl(method=trainType,
            repeats = trNumber,
            number = trBoot,
            classProbs=T,
            savePredictions = T,
            allowParallel = T)
      #print(n)
      
      result <- tryCatch({
        print (paste('train Set size =',nrow(smpl)))
        if (mdl %in% c("gbm")) {
          mFit <- train(reformulate(response=responseVar,termlabels = '.'), data=smpl,method=mdl,metric="Kappa",verbose=F, trControl=tC)
        } else {
          mFit <- train(reformulate(response=responseVar,termlabels = '.'), data=smpl,method=mdl,metric="Kappa", trControl=tC)
        }
        
        print (paste('test Set size =',nrow(testSet)))
        
        # ROC (test)
        predP <- predict(mFit,testSet,type="prob")
        cauc <- roc(testSet[[responseVar]],predP[[posClass]],auc=T)$auc+0.0
        # confusion (test)
        pred <- predict(mFit,newdata=testSet)
        conf <- confusionMatrix(pred,testSet[[responseVar]])
        cacc <- conf$overall[['Accuracy']]
        ckappa <- conf$overall[['Kappa']]
        csens <- conf$byClass[['Sensitivity']]
        cspec <- conf$byClass[['Specificity']]
        cppv <- conf$byClass[['Pos Pred Value']]
        cnpv <- conf$byClass[['Neg Pred Value']]
        cf1 <- conf$byClass[['Balanced Accuracy']]
        
        
        # ROC (train)
        predtP <- predict(mFit,newdata=smpl,type="prob")
        cauct <- roc(smpl[[responseVar]],predtP[[posClass]],auc=T)$auc
        # confusion (train)
        predt <- predict(mFit,newdata=smpl)
        conft <- confusionMatrix(predt,smpl$Diagnosis)
        cacct <- conft$overall[['Accuracy']]
        ckappat <- conft$overall[['Kappa']]
        csenst <- conft$byClass[['Sensitivity']]
        cspect <- conft$byClass[['Specificity']]
        cppvt <- conft$byClass[['Pos Pred Value']]
        cnpvt <- conft$byClass[['Neg Pred Value']]
        cf1t <- conft$byClass[['Balanced Accuracy']]        
        
        #print (paste("TRAIN: SPEC",cspec))
        #print (paste("TRAIN: SENS",csens))       
        #print (paste("TEST: SPEC",cspect))
        #print (paste("TEST: SENS",csenst))       
        
        print ('Success: training OK!')
        #print(paste("train: acc:",cacct,"kappa",ckappat,"test: acc:",cacc,'kappa:',ckappa))        
        #return(prepLearningCurveTrainTest(smpl,mdl,testSet))
        c(cacc,ckappa,cacct,ckappat,csens,cspec,cf1,cppv,cnpv,csenst,cspect,cf1t,cppvt,cnpvt,cauc,cauct)
      }, error = function(err) {
        print (paste('Error: training failed:',err))
        return(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
      })
      #print(result)
      acc <- c(acc,result[1])
      kappa <- c(kappa,result[2])
      acctrain <- c(acctrain,result[3])
      kappatrain <-c(kappatrain,result[4])
      sens <- c(sens,result[5])
      spec <- c(spec,result[6])
      f1 <- c(f1,result[7])
      ppv <- c(ppv,result[8])
      npv <- c(npv,result[9])
      senstrain <- c(senstrain,result[10])
      spectrain <- c(spectrain,result[11])
      f1train <- c(f1train,result[12])
      ppvtrain <- c(ppvtrain,result[13])
      npvtrain <- c(npvtrain,result[14])
      auc <- c(auc,result[15])
      auctrain <- c(auctrain,result[16])
      
      #print (acc)
      #print (acctrain)
      #print (kappa)
      #print (kappatrain)
      #print (paste(" -> acc:",acc,";kappa:",kappa))
      #results <- rbind(results,c(nrow(smpl),acc,kappa))
    }
    accM = mean(acc)
    accSD = sd(acc)
    kappaM = mean(kappa)
    kappaSD = sd(kappa)
    accTrainM = mean(acctrain)
    accTrainSD = sd(acctrain)
    kappaTrainM = mean(kappatrain)
    kappaTrainSD = sd(kappatrain)
    
    sensM <- mean(sens)
    sensSD <- sd(sens)
    specM <- mean(spec)
    specSD <- sd(spec)
    f1M <- mean(f1)
    f1SD <- sd(f1)
    ppvM <- mean(ppv)
    ppvSD <- sd(ppv)
    npvM <- mean(npv)
    npvSD <- sd(npv)
    sensTrainM <- mean(senstrain)
    sensTrainSD <- sd(senstrain)
    specTrainM <- mean(spectrain)
    specTrainSD <- sd(spectrain)
    f1TrainM <- mean(f1train)
    f1TrainSD <- sd(f1train)
    ppvTrainM <- mean(ppvtrain)
    ppvTrainSD <- sd(ppvtrain)
    npvTrainM <- mean(npvtrain)
    npvTrainSD <- sd(npvtrain)
    
    aucTrainM <- mean(auctrain)
    aucTrainSD <- sd(auctrain)
    aucTestM <- mean(auc)
    aucTestSD <- sd(auc)
    
    #print(paste(accM,accSD,kappaM,kappaSD))
    results <- rbind.data.frame(results,c(nrow(smpl),"Train",accTrainM,accTrainSD,kappaTrainM,kappaTrainSD,
                               sensTrainM,sensTrainSD,specTrainM,specTrainSD,f1TrainM,f1TrainSD,
                               ppvTrainM,ppvTrainSD,npvTrainM,npvTrainSD,aucTrainM,aucTrainSD))
    results <- rbind.data.frame(results,c(nrow(smpl),"Test",accM,accSD,kappaM,kappaSD,
                               sensM,sensSD,specM,specSD,f1M,f1SD,
                               ppvM,ppvSD,npvM,npvSD,aucTestM,aucTestSD))

    #print (paste(" -> Test ",cnt,'done;','Mean Acc:',accM,'; Mean Kappa:',kappaM))
    colnames(results) <- c("N","Dataset","Acc","Acc.SD","Kappa","Kappa.SD","SENS","SENS.SD","SPEC","SPEC.SD","BACC","BACC.SD","PPV","PPV.SD","NPV","NPV.SD","AUC","AUC.SD")
    #print(results)
    for (i in c(1,3:ncol(results))) {
      #print(i)
      #print(results[[i]])
      results[[i]] = as.numeric(as.character(results[[i]]))
    }
    results[[2]] <- as.character(results[[2]])    
    colnames(results) <- c("N","Dataset","Acc","Acc.SD","Kappa","Kappa.SD","SENS","SENS.SD","SPEC","SPEC.SD","BACC","BACC.SD","PPV","PPV.SD","NPV","NPV.SD","AUC","AUC.SD")
  }
  results <- as.data.frame(results)
  for (i in c(1,3:ncol(results))) {
    results[[i]] <- as.numeric(as.character(results[[i]]))
  }
  colnames(results) <- c("N","Dataset","Acc","Acc.SD","Kappa","Kappa.SD","SENS","SENS.SD","SPEC","SPEC.SD","BACC","BACC.SD","PPV","PPV.SD","NPV","NPV.SD","AUC","AUC.SD")
  results[is.na(results)] <- 0.0
  #results[is.nan(results)] <- 0.0
  results
}

plotLearningCurves <- function(lVector,tit,metrics=c("ACC","B.ACC","Kappa","NPV","PPV","Sensitivity","Specificity")) {
  vec1 <- lVector[,c(1,2,3,4)]
  vec1$Metric <- "ACC"
  colnames(vec1) <- c("N","Dataset","Value","SD","Metric")
  vec2 <- lVector[,c(1,2,5,6)]
  vec2$Metric <- "Kappa"
  colnames(vec2) <- c("N","Dataset","Value","SD","Metric")
  vec3 <- lVector[,c(1,2,7,8)]
  vec3$Metric <- "Sensitivity"
  colnames(vec3) <- c("N","Dataset","Value","SD","Metric")
  vec4 <- lVector[,c(1,2,9,10)]
  vec4$Metric <- "Specificity"
  colnames(vec4) <- c("N","Dataset","Value","SD","Metric")
  vec5 <- lVector[,c(1,2,11,12)]
  vec5$Metric <- "B.ACC"
  colnames(vec5) <- c("N","Dataset","Value","SD","Metric")
  vec6 <- lVector[,c(1,2,13,14)]
  vec6$Metric <- "PPV"
  colnames(vec6) <- c("N","Dataset","Value","SD","Metric")
  vec7 <- lVector[,c(1,2,15,16)]
  vec7$Metric <- "NPV"
  colnames(vec7) <- c("N","Dataset","Value","SD","Metric")
  vec8 <- lVector[,c(1,2,17,18)]
  vec8$Metric <- "AUC"
  colnames(vec8) <- c("N","Dataset","Value","SD","Metric")
  
  vecLong <- rbind.data.frame(vec1,vec2)
  vecLong <- rbind.data.frame(vecLong,vec3)
  vecLong <- rbind.data.frame(vecLong,vec4)
  vecLong <- rbind.data.frame(vecLong,vec5)
  vecLong <- rbind.data.frame(vecLong,vec6) #NPV
  vecLong <- rbind.data.frame(vecLong,vec7) #PPV
  vecLong <- rbind.data.frame(vecLong,vec8) #AUC
  vecLong$Metric <- as.factor(vecLong$Metric)
  vecLong$N <- as.numeric(vecLong$N)
  vecLong$Value <- as.numeric(vecLong$Value)
  vecLong$SD <- as.numeric(vecLong$SD)
  vecLong <- vecLong[vecLong$Metric %in% metrics,]
  
  #pd <- position_dodge(0.25)
  ggplot(data=vecLong,aes(x=N,y=Value,col=Dataset)) + 
    geom_line(linetype="dashed") +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Value-SD, ymax=Value+SD), colour="black", width=3) + 
    ylim(0,1) +
    xlab("Number of cases in Training set") +
    ylab("Prediction success Metric") +
    facet_grid(Metric ~ .) +
    ggtitle(tit)
}

#' Calculates prediction metrics for one model and test set
#' 
#' \code{calcPredictionMetrics} returns vector of prediction metrics
#' includes Accuracy, Sensitivity, Specificity, F1 and Kappa
#' 
#' @param testSet : test set (must be compatible with training set)
#' @param modelNames : vector of strings, names of used models, must match modles
#' @param ... : fitted models (produced by caret train)
calcPredictionMetrics <- function(mdl,testSet)
{
  pred <- predict(mdl,testSet)
  conf <- confusionMatrix(pred,testSet$Diagnosis)
  acc <- conf$overall[["Accuracy"]]
  sens <- conf$byClass[["Sensitivity"]]
  spec <- conf$byClass[["Specificity"]]
  f1 <- conf$byClass[["F1"]]
  kappa <- conf$overall[["Kappa"]]
  res <- c("Acc"=round(acc,3) ,"Kappa"=round(kappa,3),"Sens"=round(sens,3),"Spec"=round(spec,3),"F1"=round(f1,3))
  res
}

#' Calculates prediction metrics for set of models and test set
#' 
#' \code{calcPredictionMetricsTable} returns data frame of calculated ML prediction metrics
#' 
#' @param testSet : test set (must be compatible with training set)
#' @param modelNames : vector of strings, names of used models, must match modles
#' @param ... : fitted models (produced by caret train)
calcPredictionMetricsTable <- function(testSet,modelNames,...)
{
  #print (list(...))
  c = 0
  for (m in list(...)) {
    c = c + 1
    r <- calcPredictionMetrics(m,testSet)
    r <- c("Model"=modelNames[c],r)
    if (c == 1) {
      res <- r
    } else {
      res <- rbind(res,r)
    }
  }
  res <- as.data.frame(res)
  rownames(res) <- c(1:nrow(res))
  res
}

# MULTI CLASS PREDICTION FUNCTIONS
# ===========================================

#' Calculates prediction metrics and ROC curves for fitted ML/Caret multi-class prediction model, 
#' 
#' \code{generateMultiClassROC} returns list of ROC gg2plot, metrics gg2plot and metrics dataframe
#' 
#' @param fittedML : fitted Caret ML model, must enable multi-class prediction and be able to return prediction likelihood value
#' @param testSet : test set (must be compatible with training set)
#' @param tit : name of model, used for chart titles
#' @param respName : COLUMN name of response variable (to be predicted) in test set
#' @param multiROC : if TRUE, will generate averaged/multiROC curve(s)
#' @param oneAllROC : if TRUE, will generate one-vs-all ROC curves
generateMultiClassROC <- function(fittedML,testSet,tit,respName,roc.conf=T,multiROC=T,oneAllROC = T, 
                                  multiROCPriority = F, smooth=F,dottedNonPriority=T,confBoots=250,
                                  mcROCmetric = T,mcKappa = T, mcAcc = T,target = T,tarAcc = T) {
  
  # predict on test set
  pred <- predict(fittedML,testSet,type="prob")
  predAg <- as.numeric(predict(fittedML,testSet,type="raw"))
  mcAUC <- auc(multiclass.roc(testSet[[respName]], predAg))  
  
  # MAKE ROC CURVE(S):
  rocs <- list()
  rocsDFprecalc = NULL
  for (i in seq(1,ncol(pred))) {
    rocs[[i]] <- roc(testSet[[respName]]==colnames(pred)[i],pred[[i]],percent=F)
  }
  names(rocs) <- colnames(pred)
  for (i in seq(1,ncol(pred))) {
    if (class(rocsDFprecalc) != "data.frame") {
      rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
      rocsDFprecalc$outcome <- names(rocs)[i]
      rocsDFprecalc$sz = 1
    } else {
      rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
      rocsDFprecalctmp$outcome <- names(rocs)[i]
      rocsDFprecalctmp$sz = 1
      rocsDFprecalc <- rbind(rocsDFprecalc, rocsDFprecalctmp)
    }
  }
  
  # make data frame with results
  # ==========================================
  # conf matrix
  conf <- confusionMatrix(testSet[[respName]],predict(fittedML,testSet,type="raw"))
  # correct names of classes
  rownames(conf$byClass) <- gsub("Class: ","",rownames(conf$byClass))
  #rownames <- names(rocs)
  # prep results DF
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  # extract overall statistics
  predDF <- rbind.data.frame(predDF,
                             data.frame("DataSet"="multi.class","Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
  predDF <- rbind.data.frame(predDF,
                             data.frame("DataSet"="multi.class","Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
  predDF <- rbind.data.frame(predDF,
                             data.frame("DataSet"="multi.class","Metric"="ROC.AUC","Value"=round(mcAUC,2),stringsAsFactors=F))
  
  # extract per class statistics
  for (i in seq(1,nrow(conf$byClass))) {
    cc <- conf$byClass[i,]
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="Sensitivity","Value"=cc[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="Specificity","Value"=cc[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="PPV","Value"=cc[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="NPV","Value"=cc[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="B.Acc","Value"=cc[["Balanced Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="F1","Value"=cc[["F1"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="ROC.AUC","Value"=round(auc(rocs[[rownames(conf$byClass)[i]]]),2),stringsAsFactors=F))
  }
  dataSetNames <- rownames(conf$byClass)
  unik <- length(conf$byClass)
  gg2 <- ggplot(data=predDF,aes(col=Metric,y=Value,x=DataSet,shape=Metric)) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + 
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"Multi-class prediction metrics") ) +
    geom_vline(xintercept =  c(1:length(dataSetNames)) ,size=max(35,130-unik*20),alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  
  # aggregate curve to calculate mean curve
  # note: we still use non-mean curves for plot of non-mean curves
  rocCoords <- F
  rocDF <- F
  rocStep = 0.01
  for (i in seq(1,ncol(pred))) {
    if (smooth) {
      rtmp <- as.data.frame(t(coords(smooth(rocs[[i]]),x=seq(0,1,rocStep),input="sensitivity", ret = c("specificity","sensitivity") )))
    } else {
      rtmp <- as.data.frame(t(coords((rocs[[i]]),x=seq(0,1,rocStep),input="sensitivity", ret = c("specificity","sensitivity") )))
    }
    if (class(rocCoords) != "data.frame") {
      rocDF <- rtmp
      rocCoords <- rtmp
      rocDF$outcome <- names(rocs)[i]
    } else {
      rocCoords <- cbind.data.frame(rocCoords,rtmp)
      rocDFtmp <- rtmp
      rocDFtmp$outcome <- names(rocs)[i]
      rocDF <- rbind.data.frame(rocDF,rocDFtmp)
    }
  }
  rocDF$sz = 1
  rocMeanSens <- apply(rocCoords[,grep("sensitivity",colnames(rocCoords))],MARGIN = 1,mean)
  rocMeanSpec <- apply(rocCoords[,grep("specificity",colnames(rocCoords))],MARGIN = 1,mean)
  rocMDF <- data.frame(specificity=rocMeanSpec,sensitivity=rocMeanSens,outcome="mean")
  rocMDF$outcome <- as.character(rocMDF$outcome)
  rocMDF$sz = 1.1
  
  # prep dataframe for plotting
  rocDFAll <- rocsDFprecalc
  if (smooth) {
    rocDFAll <- rocDF
  }
  if (multiROC & oneAllROC) {
    rocDFAll <- rbind.data.frame(rocDFAll,rocMDF)
  } else if (multiROC & !oneAllROC) {
    rocDFAll <- rocMDF
  } else if (!multiROC & !oneAllROC) {
    print ("multiROC and oneAllROC must not all be FALSE!")
    stop()
  }
  # not doing confidence intervals
  if (!roc.conf) {
    rocDFAll <- rocDFAll[order(rocDFAll$sensitivity),]
    rocDFPlot <- rocDFAll
    rocDFPlot$ltype = rocDFPlot$sz
    if (oneAllROC & multiROC & multiROCPriority) {
      rocDFPlot$ltype[rocDFPlot$sz == 1] <- 1.25
      rocDFPlot$ltype[rocDFPlot$sz == 1.25] <- 1
    } 
    if (!dottedNonPriority) {
      rocDFPlot$ltype <- 1
    }
    gg <- ggplot(data=rocDFAll,aes(x=sensitivity,y=specificity,col=outcome)) +
      scale_size_continuous(range = c(1.25, 2)) + 
      style_roc() +
      scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
      scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
      geom_abline(slope = 1,intercept = 1,col="grey") +
      geom_line(aes(size=as.numeric(sz),linetype=as.factor(sz))) +
      guides(size = FALSE,linetype=FALSE) + #theme(legend.position="bottom") +
      ggtitle(paste(tit,"ROC") ) + xlab("Specificity") + ylab("Sensitivity") +
      labs(color='Dataset') + guides(colour = guide_legend(override.aes = list(size=3)))
  }
  # doing confidence intervals
  if (roc.conf) {
    rocDFConf <- F
    rocStep = 0.01
    boots = confBoots
    meanConf <- F
    for (i in seq(1,ncol(pred))) {
      if (smooth) {
      sens.ci <- ci.sp(smooth(rocs[[i]]), sensitivities = c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep)),
                       conf.level=0.66,boot.n = boots)
      } else {
        sens.ci <- ci.sp((rocs[[i]]), sensitivities=c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep)),
                         conf.level=0.66,boot.n = boots)
      }
      rocConf <- as.data.frame(sens.ci)
      rocConf$sensitivity <- as.numeric(rownames(rocConf))
      colnames(rocConf) <- c("sp.low","specificity","sp.high","sensitivity")
      rocConf$outcome <- names(rocs)[i]
      rocConf$sz = 1.0
      rownames(rocConf) <- NULL
      if (class(meanConf) != "data.frame") {
        meanConf <- rocConf[,c(4,1,2,3)]
      } else {
        meanConf <- cbind(meanConf,rocConf[,c(1,2,3)])
      }

      if (class(rocDFConf) != "data.frame") {
        rocDFConf <- rocConf
      } else {
        rocDFConf <- rbind.data.frame(rocDFConf,rocConf)
      }
    }
    meanConfSpLow <- apply(meanConf[,grep('sp.low',colnames(meanConf))],MARGIN = 1,mean)
    meanConfSpHigh <- apply(meanConf[,grep('sp.high',colnames(meanConf))],MARGIN = 1,mean)
    meanConfMedian <- apply(meanConf[,grep('specificity',colnames(meanConf))],MARGIN = 1,mean)
    # prep for plotting
    if (oneAllROC & multiROC & !multiROCPriority) {
      rocMDF$sp.low = 0.0
      rocMDF$sp.high = 0.0
      rocDFPlot <- rbind.data.frame(rocDFConf,rocMDF)
    } else if (oneAllROC & multiROC & multiROCPriority) {
      rocDFPlot <- rocDFConf
      rocDFPlot$sp.low = 0.0
      rocDFPlot$sp.high = 0.0
      rocDFPlotTmp <- cbind.data.frame(meanConfSpLow,meanConfMedian,meanConfSpHigh,rocConf$sensitivity)
      colnames(rocDFPlotTmp) <- c("sp.low","specificity","sp.high","sensitivity")
      rocDFPlotTmp$outcome = "multi.class"
      rocDFPlotTmp$sz = 1.25
      rocDFPlot <- rbind(rocDFPlotTmp,rocDFPlot)
    } else if (oneAllROC & !multiROC) {
      rocDFPlot <- rocDFConf
    } else if (!oneAllROC & multiROC) {
        rocDFPlot <- cbind.data.frame(meanConfSpLow,meanConfMedian,meanConfSpHigh,rocConf$sensitivity)
        colnames(rocDFPlot) <- c("sp.low","specificity","sp.high","sensitivity")
        rocDFPlot$outcome = "multi.class"
        rocDFPlot$sz = 1.25
    }
    rocDFPlot$ltype = rocDFPlot$sz
    if (oneAllROC & multiROC & multiROCPriority) {
      rocDFPlot$ltype[rocDFPlot$sz == 1] <- 1.25
      rocDFPlot$ltype[rocDFPlot$sz == 1.25] <- 1
    } 
    if (!dottedNonPriority) {
      rocDFPlot$ltype <- 1
    }
    # add 0s
    for (oc in unique(rocDFPlot$outcome)) {
      szt = rocDFPlot[rocDFPlot$outcome==oc,]$sz[1]
      ltt = rocDFPlot[rocDFPlot$outcome==oc,]$ltype[1]
      rocDFPlot <- rbind.data.frame(rocDFPlot,
                                    data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,outcome=oc,sz=szt,ltype=ltt) )
    }
    rocDFPlot <- rocDFPlot[order(rocDFPlot$specificity),]
    gg <- ggplot(rocDFPlot,aes(x=sensitivity,y=specificity,col=outcome)) +
      geom_ribbon(data=rocDFPlot,aes(ymin=sp.low,ymax=sp.high,fill=outcome),alpha=0.2,colour=NA) + 
      scale_size_continuous(range = c(1.25, 2)) + 
      style_roc() + geom_abline(slope = 1,intercept = 1,col="grey") + 
      scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
      scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
      guides(size = FALSE,linetype=FALSE) + #theme(legend.position="bottom") +
      ggtitle(paste(tit,"ROC") ) + xlab("Specificity") + ylab("Sensitivity") +
      labs(color='Dataset') + guides (fill = FALSE) + guides(colour = guide_legend(override.aes = list(size=3)))
    gg <- gg + geom_line(aes(size=as.numeric(sz),linetype=as.factor(ltype)))
  }
  
  # annotation
  # ====================
  # targets
  tarDF <- F  
  if (target) {
    for (u in unique(predDF$DataSet )) {
      if (u != "multi.class") {
        if (class(tarDF) != "data.frame") {
          tarDF <- data.frame(DataSet=u,
                              Sensitivity=predDF[predDF$DataSet==u & predDF$Metric=="Sensitivity",]$Value,
                              Specificity=predDF[predDF$DataSet==u & predDF$Metric=="Specificity",]$Value,
                              F1=predDF[predDF$DataSet==u & predDF$Metric=="F1",]$Value,
                              B.Acc=predDF[predDF$DataSet==u & predDF$Metric=="B.Acc",]$Value
                              )
        } else {
          tarDF <- rbind.data.frame(tarDF,data.frame(DataSet=u,
                                                     Sensitivity=predDF[predDF$DataSet==u & predDF$Metric=="Sensitivity",]$Value,
                                                     Specificity=predDF[predDF$DataSet==u & predDF$Metric=="Specificity",]$Value,
                                                     F1=predDF[predDF$DataSet==u & predDF$Metric=="F1",]$Value,
                                                     B.Acc=predDF[predDF$DataSet==u & predDF$Metric=="B.Acc",]$Value))
        }
      }
    }
    gg <- gg + geom_point(data=tarDF,aes(x=Specificity,y=Sensitivity,col=DataSet),size=4,shape=4,stroke=2)
  }
  
  # separate AUCses
  mvDown = 0
  if (mcROCmetric) {
    # multi-class AUC
    gg <- gg + annotate("text", x = .5, y = .11, label = paste("multi.class AUC = ",round(mcAUC,2),sep=""),hjust = 0)
    mvDown = mvDown + 1
  }
  if (mcAcc) {
    gg <- gg + annotate("text", x = .5, y = .11-0.07*mvDown, label = paste("multi.class Acc = ",round(conf$overall[["Accuracy"]],2),sep=""),hjust = 0)
    mvDown = mvDown + 1
  }
  if (mcKappa) {
    gg <- gg + annotate("text", x = .5, y = .11-0.07*mvDown, label = paste("multi.class Kappa = ",round(conf$overall[["Kappa"]],2),sep=""),hjust = 0)
  }
  if (oneAllROC) {
    posit = 0
    for (i in order(names(rocs),decreasing = T)) {
      posit = posit + 1
      cauc <- ci.auc(rocs[[i]],conf.level=0.66)
      caucsd <-  round(cauc[2]-cauc[1],2)
      if (roc.conf) {
        #lbl = paste(names(rocs)[[i]]," AUC = ",round(auc(rocs[[i]]),2),", sd = ",caucsd,sep="")
        lbl = paste(names(rocs)[[i]],": AUC = ",round(auc(rocs[[i]]),2),sep="")
      } else {
        lbl = paste(names(rocs)[[i]],": AUC = ",round(auc(rocs[[i]]),2),sep="")
      }
      if (target) {
        lbl = paste(lbl,
                    "; B.Acc = ",round(tarDF[tarDF$DataSet==names(rocs)[i], ]$B.Acc,2),
                    "; F1 = ",round(tarDF[tarDF$DataSet==names(rocs)[i], ]$F1,2),sep="")
      }
      gg <- gg + annotate("text", x = .5, y = .11+0.07*posit, label = lbl,hjust = 0)
    }
  }
  
  list(gg,gg2,predDF)
}


# normalisation function (asinsqrt, for use on metagenomes and/or pathways)
asinSqrtNormalise <- function(mN,norTaxa=T,norPWY=T) {
  rowsToNor <- c()
  if (norTaxa) {rowsToNor <- c(rowsToNor,grep('__',colnames(mN)))}
  if (norPWY) {rowsToNor <- c(rowsToNor,grep('PWY',colnames(mN)))}
  mN[,rowsToNor] <- asin(sqrt(mN[,rowsToNor]/100.0))
  mN
}
# normalisation function (divide by max & center, for use on metagenomes and/or pathways)
divideByMaxCenterNormalise <- function(mN,norTaxa=T,norPWY=T) {
  rowsToNor <- c()
  if (norTaxa) {rowsToNor <- c(rowsToNor,grep('__',colnames(mN)))}
  if (norPWY) {rowsToNor <- c(rowsToNor,grep('PWY',colnames(mN)))}
  for (r in rowsToNor) {
    mN[[r]] <- mN[[r]]-mean(mN[[r]])
    mN[[r]] <- mN[[r]]/max(abs(mN[[r]]))
  }
  mN
}
# normalisation function (divide by stdev & center, for use on metagenomes and/or pathways)
divideBySdevCenterNormalise <- function(mN,norTaxa=T,norPWY=T) {
  rowsToNor <- c()
  if (norTaxa) {rowsToNor <- c(rowsToNor,grep('__',colnames(mN)))}
  if (norPWY) {rowsToNor <- c(rowsToNor,grep('PWY',colnames(mN)))}
  for (r in rowsToNor) {
    mN[[r]] <- mN[[r]]-mean(mN[[r]])
    mN[[r]] <- mN[[r]]/sd(abs(mN[[r]]))
  }
  mN
}

# do linear correction for phenotypes
linearCorrectMGPwy <- function (iData,corrNames,corrMG=TRUE,corrPWY=TRUE,correctZeros=F,removeCorCol=T) {
  # determine which columns to correct
  #print(sum(iData==0))
  tt <- iData
  corrCols <- c()
  if (corrMG) {
    corrCols <- c(corrCols,grep("__",colnames(iData)))
  }
  if (corrPWY) {
    corrCols <- c(corrCols,grep("PWY",colnames(iData)))
  }
  # iterate over columns, correct them for selected variables (entered as corrNames)
  for (c in corrCols) {
    # build formula
    #print(c)
    if (!correctZeros) {
      corrRows <- rownames(iData)[iData[[c]] != 0.0]
    } else {
      corrRows <- rownames(iData)
    }
    toCor <- iData[corrRows,c(colnames(iData)[c], corrNames)]
    # check for factors having only 1 level:
    corCorrNames = c()
    for (cn in corrNames) {
      if (is.factor(toCor[[cn]])) {
        #print (cn)
        doDrop = F
        zer = c()
        for (lvl in levels(toCor[[cn]])) {
          if (sum(toCor[[cn]]==lvl) == 0) {
            zer = c(zer,F)
          } else {zer = c(zer,T)}
        }
        #print (zer)
        if (sum(zer) >= 2) {
          corCorrNames = c(corCorrNames,cn)
        } else {
          print(paste('WARNING:',colnames(iData)[c],': variable',cn,'has less then 2 levels! Dropping it from linear model!'))
        }
      } else {
        corCorrNames <- c(corCorrNames,cn)
      }
    }
  }
  print (corCorrNames)
  if (is.null(corCorrNames)) {
    print(paste('WARNING:',colnames(iData)[c],': Will not be corrected - all variables dropped from model'))
  } else {
      frm <- reformulate(termlabels = corCorrNames,response=colnames(iData)[c])
      # do linear model
      m <- lm(data = toCor,frm)   
      #print(m)
      # correct for linear model
      corTo <- m$coefficients[["(Intercept)"]] + resid(m)
      iData[corrRows,colnames(iData)[c]] <- corTo
  }
  # get rid of columns used in model
  if (removeCorCol) {
    for (c in corrNames) {
      iData[[c]] <- NULL
    } 
  }
  iData
}

# ===============================================
# prepares dataset by refactoring as necessary
# and taking only input covariates
# ===============================================
prepData <- function(inData,covars) {
  # grab covariates
  mdl <- inData[,covars]
  # omit NAs
  mdl <- na.omit(mdl)
  # refactor factors
  fak <- sapply(mdl, is.factor)
  c = 0
  for (i in fak) {
    c = c + 1
    if (i) {
      mdl[,c] <- as.factor(as.character(mdl[,c]))
    }
  }
  fak <- sapply(mdl, is.character)
  c = 0
  for (i in fak) {
    c = c + 1
    if (i) {
      mdl[,c] <- as.factor(as.character(mdl[,c]))
    }
  }
  mdl
}

# >>> MODEL PREP
prepModel <- function(dFrame,mdl) {
  if (mdl=="MG_PWY") {
    dFrame <- purgeMGNames(dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )])
  } else if (mdl=="P") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("P"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="P_PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("P"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="F") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("F"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="F_PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("F"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="G") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="G_PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="S") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("S"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="S_PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("S"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="GS") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("S","G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="GS_PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("S","G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(dFrame)
  } 
  colnames(dFrame) <- make.names(colnames(dFrame),unique=T)
  dFrame
}

# optimises data model using RFE, as follows:
# responseVar <- what is to be predicted [def: Diagnosis]
# trP <- percentage in training sets
# positive <- what is considered "positive"
# rfeMethod <- which method to use for optimisation [def: glm]; note: use quick one - svmRadial is also good
# testMethod <- which method to use for testing [def: svmRadial]
# xvnr <- # of x-validations [def=5]
# xvreps <- how many times to repeat x-validation [def=1]
# parallel <- if parallel processing allowed [def=T]
# verb <- prints out some info [def=T]

# RETURNS: list of
# 1: max accuracy variables
# 2: 2% tolerance in accuracy variables
# 3: 5% tolerance in accuracy variables
# 4: test set ROC for these 3 models and all-data model
# 5: x-validation ROC for these 3 models and all-data model
# 6: RFE plot
dataModelOptimiseRFE <- function(dModel,dModelName="",responseVar="Diagnosis",trP=0.75,positive="IBD",rfeMethod="glm",
                                 testMethod="svmRadial",xvnr=5,xvreps=1,parallel=T,verb=T,szs=F,trainMethod="repeatedcv") {
  
  trC <- trainControl(method="repeatedcv",number=tBut,repeats = tRep,savePredictions = T,classProbs = T,allowParallel = parallel)
  inTrain <- createDataPartition(dModel[[responseVar]],p=trP,list=F)
  trSet <- dModel[inTrain,]
  testSet <- dModel[-inTrain,]
  
  rfeCtrl <- rfeControl(functions = caretFuncs, method = trainMethod, repeats = xvreps,
                        number=xvnr, verbose = verb,allowParallel = parallel)
  if (!szs) {
    if (ncol(trSet) <= 50) {szs=seq(1,ncol(trSet)-1,1)
    } else if (ncol(trSet) <= 105) {szs=c(seq(1,50,1),seq(50,ncol(trSet)-1,2))  
    } else if (ncol(trSet) <= 210) {szs=c(seq(1,50,1),seq(50,100,2),seq(105,ncol(trSet)-1,5))
    } else if (ncol(trSet) > 210) {szs=c(seq(1,50,1),seq(50,100,2),seq(105,200,5), seq(200,ncol(trSet)-1,10))  }
  }
  #} else if (ncol(trSet) > 210) {szs=c(seq(1,50,2),seq(55,100,5),seq(120,ncol(trSet)-1,20))  }
  
  
  print(szs)
  rfeProfile <- rfe(x=trSet[,-grep(responseVar,colnames(trSet))], y=trSet$Diagnosis, sizes=szs, rfeControl = rfeCtrl,
                    metric="Kappa",method=rfeMethod)
  # select best NR of vars within tol% margin of error, while keeping number as low as possible
  varNrMax <- pickSizeTolerance(rfeProfile$results, metric="Kappa",maximize = T,tol=0.1)
  varSelMax <- caretFuncs$selectVar(y = rfeProfile$variables,size=varNrMax)
  varNr2 <- pickSizeTolerance(rfeProfile$results, metric="Kappa",maximize = T,tol=2)
  varSel2 <- caretFuncs$selectVar(y = rfeProfile$variables,size=varNr2)
  varNr5 <- pickSizeTolerance(rfeProfile$results, metric="Kappa",maximize = T,tol=6)
  varSel5 <- caretFuncs$selectVar(y = rfeProfile$variables,size=varNr5)
  # build new model with these only, compare to original model
  trC <- trainControl(method="repeatedcv",number=xvnr,repeats = xvreps,savePredictions = T,classProbs = T,allowParallel = T)
  trSetMax <- trSet[,c("Diagnosis",varSelMax)]
  trSetSel2 <- trSet[,c("Diagnosis",varSel2)]
  trSetSel5 <- trSet[,c("Diagnosis",varSel5)]
  fitAll <- train(Diagnosis ~ ., trSet,trControl = trC,method = testMethod,metric="Kappa")
  fitMax <- train(Diagnosis ~ ., trSetMax,trControl = trC,method = testMethod,metric="Kappa")
  fitRfe2 <- train(Diagnosis ~ ., trSetSel2,trControl = trC,method = testMethod,metric="Kappa")
  fitRfe5 <- train(Diagnosis ~ ., trSetSel5,trControl = trC,method = testMethod,metric="Kappa")
  # compare
  #confusionMatrix(predict(object = fitAll,newdata = trSet),trSet$Diagnosis)
  #confusionMatrix(predict(object = fitRfe2,newdata = trSet),trSet$Diagnosis)
  #confusionMatrix(predict(object = fitRfe5,newdata = trSet),trSet$Diagnosis)
  ROCfull <- roc(testSet$Diagnosis,predict(fitAll,testSet,type="prob")[[2]],auc=T,percent=F)
  ROCmax <- roc(testSet$Diagnosis,predict(fitMax,testSet,type="prob")[[2]],auc=T,percent=F)
  ROCrfe2 <- roc(testSet$Diagnosis,predict(fitRfe2,testSet,type="prob")[[2]],auc=T,percent=F)
  ROCrfe5 <- roc(testSet$Diagnosis,predict(fitRfe5,testSet,type="prob")[[2]],auc=T,percent=F)
  #g <- ggroc(list(full=ROCfull,rfe2=ROCrfe2,rfe5=ROCrfe5,max=ROCmax))
  mNames <- c(paste("All:",ncol(testSet)-1,sep=''),paste("M.max:",varNrMax,sep=''),
              paste("M.2:",varNr2,sep=''),paste("M.5:",varNr5,sep=''))
  testComp <- compareModelsOnDataset(fittedMdls = list(fitAll,fitMax,fitRfe2,fitRfe5),modelNames = mNames,
                                     posClass = positive, tit = "",mtd = paste('test',testMethod),testSet = testSet)
  cvComp <- compareModelsTrainingCV(fittedMdls = list(fitAll,fitMax,fitRfe2,fitRfe5),mtd = paste('x-val',testMethod), 
                                    modelNames = mNames,posClass = positive)
  rfeplot <- ggplot(rfeProfile) + geom_vline(xintercept = varNrMax,linetype='longdash') + geom_vline(xintercept = varNr2,linetype='longdash')+
    geom_vline(xintercept = varNr5,linetype='longdash')+ggtitle(paste("RFE plot (",dModelName,' / ',rfeMethod,')',sep=""))
  
  return(list(varSelMax,varSel2,varSel5,testComp,cvComp,rfeplot,c("Vars:Max","Vars:2%","Vars:5%","plot:testset","plot:X-val","RFE-plot")))
}

# ===============================================================================
# compares model(s) performance on one or more datasets
# ===============================================================================
compareMdlsDatasets <- function(mdls,dataSets,mdNames,posClass="IBS",roc.smooth=F,tit="",specSensAnnot=T,
                                response="Diagnosis",roc.conf = T,roc.conf.boot = 10,target = F,conf.lvl=0.66) {
  
  nrModels <- length(mdls)
  nrDataSets <- length(dataSets)
  print (nrModels)
  print (nrDataSets)
  # list of ROC curves
  rocs = list()
  # data frame of predictions
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  if (nrModels == 1 & nrDataSets == 1) {
    # ===================================================================================
    # ================== 1 MODEL, 1 DATASET =====================================
    # ===================================================================================
    leg = "Model"
    print (' -> 1 model, 1 dataset')
    trainedModels <- list()
    # generate results for each dataset
    testSet <- dataSets[[1]]
    fittedModel <- mdls[[1]]
    trainedModels[[1]] <- fittedModel
    pr <- predict(fittedModel,newdata = testSet,type="prob")
    namen = mdNames[1]
    if (roc.smooth) {
      rocs[[namen]] <- smooth(roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F))
    } else {
      rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
    }
    pr2 <- predict(fittedModel,testSet)
    conf <- confusionMatrix(pr2,testSet[[response]])
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="F1","Value"=conf$byClass[["F1"]],stringsAsFactors=F))
  } else if (nrModels == 1 & nrDataSets > 1) {
    # ===================================================================================
    # ================== 1 MODEL, MULTIPLE DATASETS =====================================
    # ===================================================================================
    print (' -> 1 model, multiple datasets')
    trainedModels <- list()
    c = 0
    # list of ROC curves
    rocs = list()
    # data frame of predictions
    predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
    # generate results for each dataset
    leg = "Dataset"
    for (ds in dataSets) {
      c = c + 1
      testSet <- dataSets[[c]]
      fittedModel <- mdls[[1]]
      trainedModels[[c]] <- fittedModel
      pr <- predict(fittedModel,newdata = testSet,type="prob")
      namen = mdNames[c]
      if (roc.smooth) {
        rocs[[namen]] <- smooth(roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F))
      } else {
        rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
      }
      pr2 <- predict(fittedModel,testSet)
      conf <- confusionMatrix(pr2,testSet[[response]])
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="F1","Value"=conf$byClass[["F1"]],stringsAsFactors=F))
    }
    
  } else if (nrModels > 1 & nrDataSets == 1) {
    # ===================================================================================
    # ================== MULTIPLE MODELs, 1 DATASET =====================================
    # ===================================================================================
    leg = "Model"
    print (' -> multiple models, 1 dataset')
    trainedModels <- list()
    c = 0
    # list of ROC curves
    rocs = list()
    # data frame of predictions
    predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
    # generate results for each dataset
    for (fittedModel in mdls) {
      c = c + 1
      testSet <- dataSets[[1]]
      trainedModels[[c]] <- fittedModel
      pr <- predict(fittedModel,newdata = testSet,type="prob")
      namen = mdNames[c]
      if (roc.smooth) {
        rocs[[namen]] <- smooth(roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F))
      } else {
        rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
      }
      pr2 <- predict(fittedModel,testSet)
      conf <- confusionMatrix(pr2,testSet[[response]])
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="F1","Value"=conf$byClass[["F1"]],stringsAsFactors=F))
    }
    
  } else if (nrModels == nrDataSets & nrModels > 1) {
    print (' -> multiple dataset - model pairs')
    trainedModels <- list()
    leg = "Data.Model"
    c = 0
    # list of ROC curves
    rocs = list()
    # data frame of predictions
    predDF <- data.frame("Data.Model"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
    # generate results for each dataset
    for (fittedModel in mdls) {
      c = c + 1
      testSet <- dataSets[[c]]
      trainedModels[[c]] <- fittedModel
      pr <- predict(fittedModel,newdata = testSet,type="prob")
      namen = mdNames[c]
      if (roc.smooth) {
        rocs[[namen]] <- smooth(roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F))
      } else {
        rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
      }
      pr2 <- predict(fittedModel,testSet)
      conf <- confusionMatrix(pr2,testSet[[response]])
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="F1","Value"=conf$byClass[["F1"]],stringsAsFactors=F))
    }
  }
  
  # extract coordinates
  if (!roc.conf) {
    rocsDFprecalc = NULL
    rocStep = 0.01
    for (i in seq(1,length(rocs))) {
      if (class(rocsDFprecalc) != "data.frame") {
        if (!roc.smooth) {
          rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalc) <- NULL
        } else {
          rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x=c(seq(0,0.1,rocStep/10),seq(0.1,1,rocStep)),ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalc) <- NULL
        }
        rocsDFprecalc[[leg]] <- names(rocs)[i]
        rocsDFprecalc$sz = 1
      } else {
        if (!roc.smooth) {
          rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalctmp) <- NULL
        } else {
          rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x=c(seq(0,0.1,rocStep/10),seq(0.1,1,rocStep)),ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalctmp) <- NULL
        }
        rocsDFprecalctmp[[leg]] <- names(rocs)[i]
        rocsDFprecalctmp$sz = 1
        rocsDFprecalc <- rbind(rocsDFprecalc, rocsDFprecalctmp)
        rownames(rocsDFprecalc) <- NULL
      }
    }
  }
  
  if (roc.conf) {
    rocDFConf <- F
    rocStep = 0.01
    boots = roc.conf.boot
    meanConf <- F
    for (i in seq(1,length(rocs))) {
      if (roc.smooth) {
        sens.ci <- ci.sp(smooth(rocs[[i]]), sensitivities =c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep)),
                         conf.level=conf.lvl,boot.n = boots)
      } else {
        sens.ci <- ci.sp((rocs[[i]]), sensitivities=c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep)),
                         conf.level=conf.lvl,boot.n = boots)
      }

      
      rocConf <- as.data.frame(sens.ci)
      rocConf$sensitivity <- as.numeric(rownames(rocConf))
      colnames(rocConf) <- c("sp.low","specificity","sp.high","sensitivity")
      rocConf[[leg]] <- names(rocs)[i]
      rocConf$sz = 1.0
      rownames(rocConf) <- NULL
      if (class(meanConf) != "data.frame") {
        meanConf <- rocConf[,c(4,1,2,3)]
      } else {
        meanConf <- cbind(meanConf,rocConf[,c(1,2,3)])
      }
      if (class(rocDFConf) != "data.frame") {
        rocDFConf <- rocConf
      } else {
        rocDFConf <- rbind.data.frame(rocDFConf,rocConf)
      }
      if (leg=="Model") {
        rocDFConf <- rbind.data.frame(rocDFConf,data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,Model=names(rocs)[i],sz=1))
      } else if (leg=="Dataset") {
        rocDFConf <- rbind.data.frame(rocDFConf,data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,Dataset=names(rocs)[i],sz=1))
      } else {
        rocDFConf <- rbind.data.frame(rocDFConf,data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,Data.Model=names(rocs)[i],sz=1))
      }
    }
    rocsDFprecalc <- rocDFConf
  }
  
  if (roc.conf) {
    rocsDFprecalc <- rocsDFprecalc[order(rocsDFprecalc$specificity,decreasing = F),]
    g <- ggplot(rocsDFprecalc,aes_string(y="specificity",x="sensitivity",col=leg)) +
      geom_ribbon(data=rocsDFprecalc,aes_string(ymin="sp.low",ymax="sp.high",fill=leg),alpha=0.2,colour=NA) + geom_line(size=1.25)
  } else {
    rocsDFprecalc <- rocsDFprecalc[order(rocsDFprecalc$sensitivity,decreasing = F),]
    g <- ggplot(rocsDFprecalc,aes_string(y="specificity",x="sensitivity",col=leg)) + geom_line(size=1.25) #pROC::ggroc(rocs) + geom_line(size=1.25)
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1) +
    ggtitle(paste(tit,"ROC") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +  
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color=leg) 
  # annotate w auc
  c = 0.0
  for (i in order(mdNames,decreasing = T)) {
    c = c + 1.0
    
    if (roc.conf) {
      rocAUCconf <- ci.auc(rocs[[ mdNames[i] ]],conf.level=conf.lvl,boot.n = boots)
      rocAUCm <- rocAUCconf[[2]]
      rocAUCdelta <- rocAUCconf[[2]] - rocAUCconf[[1]]
      g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(mdNames[i],"AUC =",round(rocAUCm,2),"\u00b1",round(rocAUCdelta,2)),hjust = 0)
    } else {
      g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(mdNames[i],"AUC =",round(rocs[[mdNames[i] ]]$auc[1],2)),hjust = 0)
    }
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  # targets
  tarDF <- F  
  if (target) {
    for (u in unique(predDF[[leg]] )) {
      if (class(tarDF) != "data.frame") {
        tarDF <- data.frame(DataSet=u,
                            Sensitivity=predDF[predDF[[leg]]==u & predDF$Metric=="Sensitivity",]$Value,
                            Specificity=predDF[predDF[[leg]]==u & predDF$Metric=="Specificity",]$Value,
                            F1=predDF[predDF[[leg]]==u & predDF$Metric=="F1",]$Value,
                            B.Acc=predDF[predDF[[leg]]==u & predDF$Metric=="B.Acc",]$Value)
      }
      g <- g + geom_point(data=tarDF,aes_string(x="Specificity",y="Sensitivity",col="DataSet"),size=4,shape=4,stroke=2)
    }
  }
  
  # do metrics plot
  unik <- max(length(dataSets),length(mdls))
  gg2 <- ggplot(data=predDF,aes_string(col="Metric",y="Value",x=leg,shape="Metric")) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + 
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"Prediction metrics") ) +
    geom_vline(xintercept =  c(1:length(mdNames)) ,size=max(35,130-unik*20),alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  list(g,gg2,predDF)
}
