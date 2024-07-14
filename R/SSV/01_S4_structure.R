precision_SSV <- methods::setClass(
  "precision_SSV",
  slots = c(
    raw.train.data = "ANY",  ## for the clustering, we only consider the train data, it should be a list
    raw.test1.data = "ANY",
    raw.test2.data = "ANY",
    harmon.train.data = "list",
    harmon.test1.data = "list",
    harmon.test2.data = "list",

    classification.result = "list",
    cluster.result = "list"
  )
)

create_precision.classification <- function(traindata, testdata1, testdata2, trainlabel, testlabel1, testlabel2) {
  object <- methods::new(Class = "precision_SSV",
                         raw.train.data = list(data = traindata, label = trainlabel),
                         raw.test1.data = list(data = testdata1, label = testlabel1),
                         raw.test2.data = list(data = testdata2, label = testlabel2)
                         )
  return(object)
}

