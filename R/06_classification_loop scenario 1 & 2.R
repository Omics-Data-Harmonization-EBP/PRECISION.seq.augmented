load("unif.train.RData")
load("unif.test.RData")
load("nonunif.train.RData")
load("nonunif.test.RData")
indexessummary=data.frame()
parallel_bench <- function(k,c1=1,c2=1,d1=0,d2=0){
  sink(tempfile())  
  on.exit(sink())
  nonuniftraindata <- nonunif.train[[k]]$data
  uniftraindata <- unif.train[[k]]$data
  
  nonuniftestdata <- nonunif.test[[k]]$data
  uniftestdata <- unif.test[[k]]$data
  
  trainlabel <- nonunif.train[[k]]$label
  testlabel <- nonunif.test[[k]]$label
  
  if(d1==0){
    traindata <- biological.effects(uniftraindata,trainlabel,c1)$data
  }
  if(d1!=0){
    traindata <- handling.effects(uniftraindata,uniftraindata,nonuniftraindata,trainlabel,d1)$data
  }
  
  if(d2==0){
    testdata <- biological.effects(uniftestdata,testlabel,c2)$data
  }
  if(d2!=0){
    testdata <- handling.effects(uniftestdata,uniftestdata,nonuniftestdata,testlabel,d2)$data
  }
  sink(tempfile())  
  on.exit(sink())
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
#The coefficient "c" represents the coefficient for biological effects, 
#while "d" represents the coefficient for handling effects.

#For scenario 1, where both the training and test groups do not include handling effects, 
#both the train data and test data should be produced using the biological.effects function, 
#and "c" should be set to 1 for both.
times=400
result_array=list()
for (k in 1:times){
  suppressMessages({
  suppressWarnings({
  result_array[[k]] <- parallel_bench(k,c1=1,c2=1,d1=0,d2=0)
  })
  })
}
tmp1 <- result_array
#For scenario 2, which mimics real-world situations, it can be assumed that both the train data and test data 
#carry a normal level of handling effect, meaning they are both produced using the handling.effects function, 
#and "d" should be set to 1 for both. 
times=400
result_array=list()
for (k in 1:times){
  suppressMessages({
  suppressWarnings({
  result_array[[k]] <- parallel_bench(k,c1=1,c2=1,d1=1,d2=1)
  })
  })
}
tmp2 <- result_array
#Another assumption we discussed in the meeting is that the train data is clean, 
#without handling effects, while the test data carries a normal level of handling effects. 
#This means the train data is produced using the biological.effects function, 
#while the test data is produced using the handling.effects function.
times=400
result_array=list()
for (k in 1:times){
  suppressMessages({
  suppressWarnings({
  result_array[[k]] <- parallel_bench(k,c1=1,c2=1,d1=0,d2=1)
  })
  })
}
tmp3 <- result_array

