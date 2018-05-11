#Started on Monday May 7, 2018 by Humza Haider. The goal of this exploratory file 
#is to discern how to retrieve the survival curves from each node of each tree in 
#a random forest. 
library(randomForestSRC);
library(caret)
library(data.table)
library(zoo)
library(parallel)
library(survival)

runRSF = function(dataset, timeVar, censVar,numCV, ntree = 1000, nsplit = 3, nodesize =3){
  set.seed(42)
  names(dataset)[which(names(dataset) == timeVar)] = "time"
  names(dataset)[which(names(dataset) == censVar)] = "status"
  foldIndex = createFolds(dataset$status, k = numCV, list = F)
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
        summ = summary(z[[1]])
        cbind.data.frame(TreeID = x, NodeID = z[[2]],time = c(0,summ$time), survival = c(1,summ$surv), N = z[[1]]$n, numEvents = sum(z[[1]]$n.event))
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
    patientList = data.table(ID = as.numeric(),Time=as.numeric(),Survival=as.numeric())
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
      setkey(DTFull, TreeID, cutTime)
      #Adds the most recent survival probability to NA values going forward.
      DTFull[,totalEvents := sum(DTFull[DTFull[,.I[1], c('TreeID','numEvents')]$V1,numEvents], na.rm=T)]
      DTFull[,survival := survival*(numEvents/totalEvents)]
      
      DTFull[,survival := na.locf(survival)]
      IndividualCurve = DTFull[, sum(survival), by =cutTime]
      dats = data.table(ID = i,Time = times[-length(times)], Survival = IndividualCurve[,V1])
      patientList = rbindlist(list(patientList,dats))
    }
    patientList[,CVID := CV]
    toReturn = rbind.data.frame(toReturn,patientList)
  }
  return(toReturn)
}

smallPat = patientList[patientList$ID %in% 40:60,]
smallPat$ID = as.factor(smallPat$ID)
ggplot(smallPat,aes(Time,Survival,color=ID, group =ID))+
        geom_point()+
        geom_line()+
        theme_bw()+
        theme(text = element_text(size = 15,face = "bold"))+
        labs(y = "Survival Probability","Time")



