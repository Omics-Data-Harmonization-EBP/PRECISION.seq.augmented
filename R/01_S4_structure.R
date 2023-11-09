precision <- methods::setClass(
  "precision",
  slots = c(
    raw.train.data = "list",
    raw.test.data = "list",
    harmon.train.data = "list",
    harmon.test.data = "list",
    classification.train.result = "list",
    classification.test.result = "list",
    raw.data = "list",
    harmon.data = "list",
    cluster.result = "list"
  )
)

create_precision.cluster <- function(data, label) {
  object <- methods::new(Class = "precision",
                         raw.data = list(data = data, label = label))
  return(object)
}

create_precision.classification <- function(traindata, testdata, trainlabel, testlabel) {
  object <- methods::new(Class = "precision",
                         raw.train.data = list(data = traindata, label = trainlabel),
                         raw.test.data = list(data = testdata, label = testlabel)
                         )
  return(object)
}
