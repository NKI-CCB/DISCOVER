import numpy

from _discover import maxent


def count_unique(keys):
    # From http://stackoverflow.com/questions/10741346/numpy-frequency-counts-for-unique-values-in-an-array
    uniq_keys = numpy.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, numpy.bincount(bins)


def estimateBackground(events, strata=None):
    events = numpy.asarray(events)
   
    if strata is None:
        rowSums = events.sum(1)
        colSums = events.sum(0)

        rowValues, rowWeights = count_unique(rowSums)
        colValues, colWeights = count_unique(colSums)

        mu = maxent.fit(rowValues, rowWeights, colValues, colWeights)
        nRows = len(rowValues)
        eA = numpy.dot(numpy.exp((mu[:nRows] / rowWeights)[:, numpy.newaxis]),
                       numpy.exp((mu[nRows:] / colWeights)[numpy.newaxis]))
        P = 1.0 / (eA + 1)

        if P.max() > 1:
            import warnings
            warnings.warn("Some background estimates are greater than 1")
            #P[P > 1] = 1

        return P[rowValues.searchsorted(rowSums)[:, numpy.newaxis],
                 colValues.searchsorted(colSums)[numpy.newaxis]]
    else:
        strata = numpy.asarray(strata)
        result = numpy.empty(events.shape)
        for x in numpy.unique(strata):
            result[:, strata == x] = estimateBackground(events[:, strata == x])
        return result
