#Started on Monday May 7, 2018 by Humza Haider. The goal of this exploratory file 
#is to discern how to retrieve the survival curves from each node of each tree in 
#a random forest. 

library(randomForestSRC);
library(caret)
library(data.table)
library(zoo)
library(parallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores, type = "FORK")

data(wihs, package = "randomForestSRC")
wihs = wihs[-which(wihs$status ==2),]
inTrain = createDataPartition(wihs$time, p =.7,list=F)
train = wihs[inTrain,]
test = wihs[-inTrain,]
ntree = 100
set.seed(42)
wihs.obj <- rfsrc(Surv(time, status) ~ ., train, nsplit = 3, ntree = ntree,
                  membership = T)
#Get the training nodes
trainMembership = as.data.frame(wihs.obj$membership)
listOfCurves =  parLapply(cl,1:ntree, function(x){
                          treeInfo = trainMembership[,x]
                          uniqueNodes= unique(treeInfo)
                          listOfStuff = lapply(uniqueNodes, function(y){summ = summary(survfit(Surv(time,status)~1, train[which(treeInfo == y),]));
                          node = y;
                          list(summ, node)})
                          plyr::ldply(listOfStuff, function(z) cbind.data.frame(TreeID = x, NodeID = z[[2]],time = c(0,z[[1]]$time), survival = c(1,z[[1]]$surv)))
                          })
KMCurves = plyr::ldply(listOfCurves, rbind)

DT = data.table(KMCurves)
setkey(DT, TreeID,NodeID)

predictions = predict(wihs.obj, test, membership = T)
#Getting the tree node information:
nodeMembership = as.data.frame(predictions$membership)
predictMe = data.table(t(nodeMembership))
patientList = data.table(ID = as.numeric(),Time=as.numeric(),Survival=as.numeric())
for(i in 1:(ncol(predictMe))){
  person = predictMe[,..i]
  personCurves = DT[.(1:100, person)]
  times = personCurves[,sort(unique(time))]
  
  #We add some time at the end so the last two time points don't end up in the same interval.
  times =c(times,times[length(times)] + 1e-6)
  cuttedTimes = cut(times[-62],times,include.lowest = T,right=F)
  personCurves[, cutTime := cut(time, times, include.lowest = T,right=F)]
  
  toJoin = data.table(TreeID = rep(1:100,each = length(cuttedTimes)),cutTime = cuttedTimes)
  DTFull = merge(toJoin, personCurves, by = c("TreeID","cutTime"), all.x = T)
  
  #Sets up ordering for na.locf
  setkey(DTFull, TreeID, cutTime)
  #Adds the most recent survival probability to NA values going forward.
  DTFull[,survival := na.locf(survival)]
  IndividualCurve = DTFull[, mean(survival), by =cutTime]
  dats = data.table(ID = i,Time = times[-length(times)], Survival = IndividualCurve[,V1])
  patientList = rbindlist(list(patientList,dats))
  if(i %% round(ncol(predictMe)/10) == 0){
    print(paste(i*100/ncol(predictMe),"%",sep = ""))
  }
}
smallPat = patientList[patientList$ID %in% 1:10,]
smallPat$ID = as.factor(smallPat$ID)
ggplot(smallPat,aes(Time,Survival,color=ID, group =ID))+
        geom_point()+
        geom_line()+
        theme_bw()+
        theme(text = element_text(size = 15,face = "bold"))+
        labs(y = "Survival Probability","Time")



