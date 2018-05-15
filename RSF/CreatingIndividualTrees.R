#Started on Monday May 7, 2018 by Humza Haider. The goal of this exploratory file 
#is to discern how to retrieve the survival curves from each node of each tree in 
#a random forest. 
library(randomForestSRC);
library(caret)
library(data.table)
library(zoo)
library(parallel)
library(survival)

runAvgRSF = function(dataset, timeVar, censVar,numCV, ntree = 1000, nsplit = 3, nodesize =3){
  set.seed(42)
  names(dataset)[which(names(dataset) == timeVar)] = "time"
  names(dataset)[which(names(dataset) == censVar)] = "status"
  foldIndex = createFolds(as.factor(dataset$status), k = numCV, list = F)
  toReturn = data.frame()
  for(CV in 1:numCV){
    test = dataset[foldIndex == CV,]
    train = dataset[foldIndex!=CV,,]
    forest = rfsrc(Surv(time, status) ~ ., train, nsplit = nsplit, ntree = ntree, nodesize = nodesize, membership = T)
    trainMembership = as.data.frame(forest$membership)
    no_cores <- detectCores()
    cl <- makeCluster(no_cores-1, type = "FORK")
    listOfCurves =  parLapply(cl,1:ntree, function(x){
      treeInfo = trainMembership[,x]
      uniqueNodes= unique(treeInfo)
      listOfStuff = lapply(uniqueNodes, function(y){
        survMod = survfit(Surv(time,status)~1, train[which(treeInfo == y),])
        node = y
        list(survMod, node)})
      plyr::ldply(listOfStuff, function(z){
        mod = z[[1]]
        cbind.data.frame(TreeID = x, NodeID = z[[2]],time = c(0,mod$time), survival = c(1,mod$surv), N = c(0,mod$n.event),
                         numCens = c(0,mod$n.cens))
      })
    })
    stopCluster(cl)
    KMCurves = plyr::ldply(listOfCurves, rbind)
    
    DT = data.table(KMCurves)
    setkey(DT, TreeID,NodeID)
    
    predictions = predict(forest, test, membership = T)
    #Getting the tree node information:
    nodeMembership = as.data.frame(predictions$membership)
    predictMe = data.table(t(nodeMembership))
    patientList = data.table(ID = as.numeric(),OriginalIndex = as.numeric(),Time=as.numeric(),Survival=as.numeric())
    indexList = which(foldIndex == CV)
    for(i in 1:(ncol(predictMe))){
      person = predictMe[,..i]
      personCurves = DT[.(1:ntree, person)]
      times = personCurves[,sort(unique(time))]
      
      #We add some time at the end so the last two time points don't end up in the same interval.
      times =c(times,times[length(times)] + 1e-6)
      cuttedTimes = cut(times[-length(times)],times,include.lowest = T,right=F)
      personCurves[, cutTime := cut(time, times, include.lowest = T,right=F)]
      
      toJoin = data.table(TreeID = rep(1:ntree,each = length(cuttedTimes)),cutTime = cuttedTimes)
      DTFull = merge(toJoin, personCurves, by = c("TreeID","cutTime"), all.x = T)
      
      #Sets up ordering for na.locf
      setkey(DTFull, cutTime)
      
      #Adds the most recent survival probability to NA values going forward.
      DTFull[,N := ifelse(is.na(N),0,N)]
      DTFull[,numCens :=ifelse(is.na(numCens),0,numCens)]
      
      DTFull[,totalEvents := length(which(time!=0 & !is.na(time))),by =  c('TreeID')]
      
      DTFull[,runningCens :=c(0,cumsum(numCens[-.N])), by = "TreeID"]
      DTFull[,runningAllCens :=c(0,cumsum(numCens[-.N]))]
      
      DTFull[,runningEvents := sum(DTFull[DTFull[,.I[1], c('TreeID','totalEvents')]$V1,totalEvents], na.rm=T) - runningAllCens]
      DTFull[,totalEvents2 := totalEvents - runningCens]
      
      DTFull[,survival := ifelse(time!=0 &N!=1, NA,survival)]
      
      
      setkey(DTFull, TreeID, cutTime)
      DTFull[,survival := na.locf(survival)]
      setkey(DTFull, cutTime)
      
      DTFull[,weightInd := ifelse(sum(N) >0,1,0), by = "cutTime"]
      DTFull[,survival := ifelse(weightInd | time ==0,survival*(totalEvents2/runningEvents),NA )]
      setkey(DTFull, TreeID, cutTime)
      DTFull[,survival := na.locf(survival)]

      
      #Adds the most recent survival probability to NA values going forward.
      IndividualCurve = DTFull[, sum(survival), by =cutTime]
      dats = data.table(ID = i,OriginalIndex = indexList[i],Time = times[-length(times)], Survival = IndividualCurve[,V1])
      patientList = rbindlist(list(patientList,dats))
    }
    patientList[,CVID := CV]
    toReturn = rbind.data.frame(toReturn,patientList)
  }
  return(toReturn)
}


runCumRSF = function(dataset, timeVar, censVar,numCV, ntree = 1000, nsplit = 3, nodesize =3){
  set.seed(42)
  names(dataset)[which(names(dataset) == timeVar)] = "time"
  names(dataset)[which(names(dataset) == censVar)] = "status"
  foldIndex = createFolds(as.factor(dataset$status), k = numCV, list = F)
  toReturn = data.frame()
  for(CV in 1:numCV){
    test = dataset[foldIndex == CV,]
    train = dataset[foldIndex!=CV,,]
    forest = rfsrc(Surv(time, status) ~ ., train, nsplit = nsplit, ntree = ntree, nodesize = nodesize, membership = T)
    trainMembership = as.data.frame(forest$membership)
    no_cores <- detectCores()
    cl <- makeCluster(no_cores-1, type = "FORK")
    listOfCurves =  parLapply(cl,1:ntree, function(x){
      treeInfo = trainMembership[,x]
      uniqueNodes= unique(treeInfo)
      listOfStuff = lapply(uniqueNodes, function(y){
        survMod = survfit(Surv(time,status)~1, train[which(treeInfo == y),])
        node = y
        list(survMod, node)})
      plyr::ldply(listOfStuff, function(z){
        mod = z[[1]]
        cbind.data.frame(TreeID = x, NodeID = z[[2]],time = mod$time, Nevent = mod$n.event, Ncens = mod$n.censor)
      })
    })
    stopCluster(cl)

    KMCurves = plyr::ldply(listOfCurves, rbind)
    
    DT = data.table(KMCurves)
    setkey(DT, TreeID,NodeID)
    
    predictions = predict(forest, test, membership = T)
    #Getting the tree node information:
    nodeMembership = as.data.frame(predictions$membership)
    predictMe = data.table(t(nodeMembership))
    patientList = data.table(ID = as.numeric(),OriginalIndex = as.numeric(),Time=as.numeric(),Status=as.numeric())
    indexList = which(foldIndex == CV)
    for(i in 1:(ncol(predictMe))){
      person = predictMe[,..i]
      personCurves = DT[.(1:ntree, person)]
      
      personCurvesUnc= personCurves[Nevent>0,]
      personCurvesCen= personCurves[Ncens>0,]
      
      personCurvesUnc[,Ncens :=0]
      personCurvesCen[,Nevent :=0]
      
      personCurvesUncExpanded = personCurvesUnc[rep(seq_len(nrow(personCurvesUnc)), personCurvesUnc$Nevent),1:5]
      personCurvesCenExpanded = personCurvesCen[rep(seq_len(nrow(personCurvesCen)), personCurvesCen$Ncens),1:5]

      personCurvesExpanded = rbindlist(list(personCurvesCenExpanded, personCurvesUncExpanded))
      #Make indicator for death vs censoring. Need to catch when death and censoring occur in same time point.
      #personCurvesExpanded[,totalEvents := Nevent+Ncens, by = c("TreeID","NodeID","time")]
      personCurvesExpanded[,status := ifelse(Nevent,1,0) ,by = c("TreeID","NodeID","time")]
      
      dats = data.table(ID = i,OriginalIndex = indexList[i],Time = personCurvesExpanded[,time], Status = personCurvesExpanded[,status])
      patientList = rbindlist(list(patientList,dats))
    }
    patientList[,CVID := CV]
    toReturn = rbind.data.frame(toReturn,patientList)
  }
  return(toReturn)
}

banana = runAvgRSF(lung, "time","staus",ntree = 2,numCV = 3)
banana2 = runCumRSF(lung, "time","staus",ntree = 2,numCV = 3)

smallPat = banana[banana$ID %in% j & banana$CVID == 1]
smallPat$ID = as.factor(smallPat$ID)
ggplot(smallPat,aes(Time,Survival,color=ID, group =ID))+
        geom_point()+
        geom_line()+
        theme_bw()+
        theme(text = element_text(size = 15,face = "bold"))+
        labs(y = "Survival Probability","Time")

smallPat2 = banana2[ID == j & CVID == 1]
plot(survfit(Surv(smallPat2$Time, smallPat2$Status)~1))

smallPat$Survival
c(1,summary(survfit(Surv(smallPat2$Time, smallPat2$Status)~1))$surv)
smallPat2
j=j+1



library(ipred)

getIndividualBrier = function(testList, individualID, predictionCurves){
  curve = predictionCurves[ID == individualID,]
  
}

testing = unique(banana[CVID ==1,OriginalIndex])
orig = lung[testing,]
smod = Surv(orig$time, orig$status)
sbrier()













