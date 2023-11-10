precision <- methods::setClass(
  "precision",
  slots = c(
    raw.train.data = "list",  ## for the clustering, we only consider the train data
    raw.test.data = "list",
    harmon.train.data = "list",
    harmon.test.data = "list",

    classification.result = "list",
    cluster.result = "list"
  )
)

create_precision.cluster <- function(data, label) {
  object <- methods::new(Class = "precision",
                         raw.train.data = list(data = data, label = label),
                         raw.test.data = NULL)
  return(object)
}

create_precision.classification <- function(traindata, testdata, trainlabel, testlabel) {
  object <- methods::new(Class = "precision",
                         raw.train.data = list(data = traindata, label = trainlabel),
                         raw.test.data = list(data = testdata, label = testlabel)
                         )
  return(object)
}
