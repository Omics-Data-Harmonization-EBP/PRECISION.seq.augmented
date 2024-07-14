#let's say we have generated three dataset: train/test1/test2
load("train.dataset.RData")
load("test1.dataset.RData")
load("test2.dataset.RData")

indexessummary=data.frame()
parallel_bench <- function(k){
  sink(tempfile())  
  on.exit(sink())
  traindata <- train.dataset[[k]]$data
  testdata1 <- test1.dataset[[k]]$data
  testdata2 <- test2.dataset[[k]]$data
  
  trainlabel <- train.dataset[[k]]$label
  testlabel1 <- test1.dataset[[k]]$label
  testlabel2 <- test2.dataset[[k]]$label
  
  object1 <- create_precision.classification(traindata=traindata,testdata1=testdata1,testdata2=testdata2,
                                             trainlabel=trainlabel,testlabel1=testlabel1,testlabel2=testlabel2)
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
  return(object1@classification.result)

}

times=1
result_array=list()
for (k in 1:times){
  suppressMessages({
  suppressWarnings({
  result_array[[k]] <- parallel_bench(k)
  })
  })
}


