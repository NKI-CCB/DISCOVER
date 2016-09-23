import numpy
import pandas

from collections import namedtuple
from _discover import fdr


def pairwise_discover_test(x, g=None, alternative="less", correct=True):
    """
    Perform many pairwise mutual exclusivity or co-occurrence tests.

    Parameters
    ----------
    x : DiscoverMatrix

    g : array_like, optional
        An optional grouping vector for the rows of `x`. Pairs of rows within
        the same group are not tested.

    alternative : {'less', 'greater'}, optional
        If 'less', a mutual-exclusivity analysis is performed, if 'greater' a
        co-occurrence analysis.

    correct : bool, optional
        If True, multiple testing correction is performed.

    Returns
    -------
    result : PairwiseDiscoverResult
        An object containing the test results for all pairwise combinations.
    """
    assert alternative in ["less", "greater"]

    events = x.events
    bg = x.bg

    if g is None:
        pFlat, qFlat, pi0 = fdr.mutex(events, bg, alternative == "less")
        
        p = numpy.empty((x.shape[0], ) * 2)
        p[:] = numpy.nan
        p[numpy.triu_indices_from(p, 1)] = pFlat

        q = numpy.empty((x.shape[0], ) * 2)
        q[:] = numpy.nan
        q[numpy.triu_indices_from(p, 1)] = qFlat
    else:
        i = numpy.argsort(g)
        levels, inverse = numpy.unique(g, return_inverse=True)
        blockSizes = numpy.bincount(inverse)
        p, q, pi0 = fdr.analyseblockstructure(events[i], bg[i], alternative == "less", blockSizes)

        j = numpy.argsort(i)
        p = p[j[:, numpy.newaxis], j]
        q = q[j[:, numpy.newaxis], j]

    p = pandas.DataFrame(p, index=x.rownames, columns=x.rownames)
    q = pandas.DataFrame(q, index=x.rownames, columns=x.rownames)

    return PairwiseDiscoverResult(p, q, pi0, alternative)


class PairwiseDiscoverResult:
    """
    Class to store the results of the `pairwise_discover_test` function.

    Attributes
    ----------
    pvalues, qvalues : pandas.DataFrame
        Matrices containing the pairwise DISCOVER test P values and
        FDR-corrected Q values. The P value corresponding to a gene pair
        (`gene1`, `gene2`) is stored in either `pvalues.ix[gene1, gene2]` or
        `pvalues.ix[gene2, gene1]`. The other entry will contain ``nan``.
        Q values are stored in the same way.

    pi0 : float
        Estimate of the proportion of true null hypotheses.

    alternative : {'less', 'greater'}
        If 'less', these results relate to mutual exclusivity, if 'greater' to
        co-occurrence.
    """

    def __init__(self, pvalues, qvalues, pi0, alternative):
        self.pvalues = pvalues
        self.qvalues = qvalues
        self.pi0 = pi0
        self.alternative = alternative

    def significant_pairs(self, q_threshold=0.01):
        """
        Return the gene pairs significant at a specified maximum false discovery
        rate.

        Parameters
        ----------
        q_threshold : float, optional
            The maximum false discovery rate (default 0.01) for test results
            included in the result.

        Returns
        -------
        result : pandas.DataFrame
            DataFrame with significant gene pairs. The gene names are stored in
            the columns ``gene1`` and ``gene2``. The columns ``pvalue`` and
            ``qvalue`` contain the DISCOVER test P value and the FDR-corrected
            Q value.
        """
        i, j = numpy.where(self.pi0 * numpy.asarray(self.qvalues) < q_threshold)
        return pandas.DataFrame({
            "gene1": self.qvalues.index[i],
            "gene2": self.qvalues.columns[j],
            "pvalue": numpy.asarray(self.pvalues)[i, j],
            "qvalue": numpy.asarray(self.qvalues)[i, j]})

    def __repr__(self):
        return (
            "Pairwise DISCOVER {type} test\n"
            "alternative hypothesis: observed overlap is {alternative} than expected by chance\n"
            "\n"
            "number of pairs tested: {num_tested}\n"
            "proportion of true null hypotheses: {pi0}\n"
            "number of significant pairs at a maximum FDR of {fdr_threshold}: {num_significant}\n").format(
                type={"less": "mutual exclusivity", "greater": "co-occurrence"}[self.alternative],
                alternative=self.alternative,
                num_tested=numpy.isfinite(numpy.asarray(self.pvalues)).sum(),
                pi0=self.pi0,
                fdr_threshold=0.01,
                num_significant=numpy.sum(self.pi0 * numpy.asarray(self.qvalues) < 0.01))
