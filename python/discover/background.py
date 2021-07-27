import numpy

from _discover import maxent


def count_unique(keys):
    # From http://stackoverflow.com/questions/10741346/numpy-frequency-counts-for-unique-values-in-an-array
    uniq_keys = numpy.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, numpy.bincount(bins)


def estimate_background(events, strata=None):
    events = numpy.asarray(events)
   
    if strata is None:
        rowSums = events.sum(1)
        colSums = events.sum(0)

        rowValues, rowWeights = count_unique(rowSums)
        colValues, colWeights = count_unique(colSums)

        mu = maxent.fit(rowValues, rowWeights, colValues, colWeights)
        nRows = len(rowValues)
        eA = numpy.outer(numpy.exp(mu[:nRows]), numpy.exp(mu[nRows:]))
        P = 1.0 / (eA + 1)

        return P[rowValues.searchsorted(rowSums)[:, numpy.newaxis],
                 colValues.searchsorted(colSums)[numpy.newaxis]]
    else:
        strata = numpy.asarray(strata)
        result = numpy.empty(events.shape)
        for x in numpy.unique(strata):
            result[:, strata == x] = estimate_background(events[:, strata == x])
        return result
