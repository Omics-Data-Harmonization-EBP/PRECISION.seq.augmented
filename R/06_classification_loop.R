load("unif.train.RData")
load("unif.test.RData")
load("nonunif.train.RData")
load("nonunif.test.RData")
indexessummary=data.frame()
parallel_bench <- function(k,c=1,d=1){
  
  nonuniftraindata <- nonunif.train[[k]]$data
  uniftraindata <- unif.train[[k]]$data
  testdata <- unif.test[[k]]$data
  trainlabel <- nonunif.train[[k]]$label
  testlabel <- nonunif.test[[k]]$label
  
  
  
  traindata <- handling.effects(uniftraindata,uniftraindata,nonuniftraindata,trainlabel,d)$data
  testdata <- biological.effects(testdata,testlabel,c)$data
  
  object1 <- create_precision.classification(traindata=traindata,testdata=testdata,trainlabel=trainlabel,testlabel=testlabel)
  object1 <- harmon.TC(object1)
  object1 <- harmon.UQ(object1)
  object1 <- harmon.med(object1)
  object1 <- harmon.TMM(object1)
  object1 <- harmon.DESeq(object1)
  object1 <- harmon.PoissonSeq(object1)
  object1 <- harmon.QN(object1)
  object1 <- harmon.RUVg(object1)
  object1 <- harmon.RUVs(object1)
  object1 <- harmon.RUVr(object1)
  object1 <- classification.pam(object1)
  object1 <- classification.knn(object1)
  object1 <- classification.lasso(object1)
  object1 <- classification.svm(object1)
  object1 <- classification.ranfor(object1)
  object1 <- classification.all(object1)
  return(object1@classification.result)
}

times=400
result_array=list()
for (k in 1:times){
  result_array[[k]] <- parallel_bench(k,c=1,d=1)
}

#length(result_array)
save(result_array,file = "classification result.RData")