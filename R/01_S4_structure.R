precision <- methods::setClass(
  "precision",
  slots = c(
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

