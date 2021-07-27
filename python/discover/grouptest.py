import numpy

from _discover import poisbinom


def groupTestCoverage(events, bg):
    pCovered = 1 - numpy.exp(numpy.log1p(-bg).sum(0))
    return 1 - poisbinom.cdf(pCovered, events.any(0).sum() - 1)


def groupTestExclusivity(events, bg):
    logP = numpy.log(bg)
    logNotP = numpy.log(1 - bg)
    pExactlyOne = numpy.exp(
        numpy.logaddexp.reduce(
            numpy.repeat(logNotP.sum(0)[numpy.newaxis], bg.shape[0], 0) + logP - logNotP, 0))
    x = numpy.sum(events.sum(0) == 1)

    return 1 - poisbinom.cdf(pExactlyOne, x - 1)


def groupTestImpurity(events, bg):
    logP = numpy.log(bg)
    logNotP = numpy.log(1 - bg)
    pNone = numpy.exp(logNotP.sum(0)[numpy.newaxis])
    pExactlyOne = numpy.exp(
        numpy.logaddexp.reduce(
            numpy.repeat(logNotP.sum(0)[numpy.newaxis], bg.shape[0], 0) + logP - logNotP, 0))

    pMoreThanOne = 1 - pNone - pExactlyOne
    # Due to floating point precision issues the above probabilities
    # might end up being slightly smaller than 0; clip them at 0
    pMoreThanOne = numpy.maximum(0, pMoreThanOne)

    x = numpy.sum(events.sum(0) > 1)
    return poisbinom.cdf(pMoreThanOne, x)


def groupwise_discover_test(events, method="impurity"):
    """
    Perform a groupwise mutual exclusivity test.

    Parameters
    ----------
    events : DiscoverMatrix
        Matrix with rows corresponding to the genes in the gene set to be
        tested.

    method : {'impurity', 'coverage', 'exclusivity'}
        The mutual exclusivity statistic to estimate significance for.

    Returns
    -------
    pvalue : float
        The P value of the groupwise DISCOVER test.
    """
    methods = {
        "coverage": groupTestCoverage,
        "exclusivity": groupTestExclusivity,
        "impurity": groupTestImpurity
    }
    
    try:
        return methods[method](events.events, events.bg)
    except KeyError:
        raise ValueError("Unknown method: %s" % method)
