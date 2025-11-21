# Precision Class

This class handles the data and results for the PRECISION.seq.augmented
package. It includes raw training and test data, harmonized data for
various normalization methods, classification results, and clustering
results. Use
[`create.precision.classification`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/create.precision.classification.md)
to create an instance for classification tasks, and
[`create.precision.cluster`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/create.precision.cluster.md)
to create an instance for clustering tasks.

## Slots

- `raw.train.data`:

  A list containing the training data (`raw.train.data$data`) and labels
  (`raw.train.data$label`)

- `raw.test1.data`:

  A list containing the first test dataset (`raw.test1.data$data`) and
  labels (`raw.test1.data$label`)

- `raw.test2.data`:

  A list containing the second test dataset (`raw.test2.data$data`) and
  labels (`raw.test2.data$label`)

- `harmon.train.data`:

  A list of harmonized training data `raw.train.data` for various
  methods

- `harmon.test1.data`:

  A list of harmonized first test data `raw.test1.data` for various
  methods

- `harmon.test2.data`:

  A list of harmonized second test data `raw.test2.data` for various
  methods

- `classification.result`:

  A list containing the results of classification tasks

- `cluster.result`:

  A list containing the results of clustering tasks

## See also

[`create.precision.classification`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/create.precision.classification.md)
for creating instances of this class for classification,
[`create.precision.cluster`](https://omics-data-harmonization-ebp.github.io/PRECISION.seq.augmented/reference/create.precision.cluster.md)
for creating instances of this class for clustering.
