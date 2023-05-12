import numpy
import pandas

from .background import estimate_background


class DiscoverMatrix(object):
    """
    Class to store a binary alteration matrix and the corresponding alteration
    probability estimates.

    Parameters
    ----------
    events : pandas.DataFrame
        Binary alteration matrix.

    bg : pandas.DataFrame, optional
        Matrix of alteration probabilities for `events`. Mainly for internal use
        only. Use at your own risk.

    strata : array_like, optional
        To perform a stratified DISCOVER test, this array should contain the
        strata. The length of this array must match the number of columns of
        `events`.

    Attributes
    ----------
    events : ndarray
        Binary alteration matrix.

    bg : ndarray
        Matrix with estimated (background) alteration probabilities.

    rownames : ndarray
        Array containing the row names for `events` and `bg`.

    colnames : ndarray
        Array containing the column names for `events and `bg`.

    shape : tuple of int
        The shape of `events` and `bg`.
    """

    def __init__(self, events, bg=None, strata=None):
        if not isinstance(events, pandas.DataFrame):
            raise TypeError("Parameter `events` must be a pandas.DataFrame")

        self._events = events.copy(deep=False)

        if bg is None:
            self._bg = pandas.DataFrame(estimate_background(events, strata), index=events.index, columns=events.columns)
        else:
            if strata is not None:
                import warnings
                warnings.warn("The value of 'strata' is ignored if 'bg' is provided.")

            self._bg = bg.copy(deep=False)
       
        self.check_consistency()

    def check_consistency(self):
        assert self._events.shape == self._bg.shape
        assert (self._events.index == self._bg.index).all()
        assert (self._events.columns == self._bg.columns).all()

    def __getitem__(self, idx):
        if not isinstance(idx, (pandas.Index, pandas.Series)):
            idx = numpy.asanyarray(idx)

        if not idx.ndim == 1:
            raise ValueError("`idx` should be a 1-dimensional array")

        # This is meant to emulate pandas' deprecated .ix[] indexer.
        # If `idx` contains integers, use .iloc[] unless `_events.index`
        # has an integer dtype. In that case, as well as when `idx` does
        # not contain integers, use .loc[].
        if (numpy.issubdtype(idx.dtype, numpy.integer) and
              not numpy.issubdtype(self._events.index.dtype, numpy.integer)):
            return DiscoverMatrix(self._events.iloc[idx], self._bg.iloc[idx])
        else:
            return DiscoverMatrix(self._events.loc[idx], self._bg.loc[idx])

    @property
    def events(self):
        return numpy.asarray(self._events)

    @property
    def bg(self):
        return numpy.asarray(self._bg)

    @property
    def rownames(self):
        return numpy.asarray(self._events.index)

    @rownames.setter
    def rownames(self, names):
        self._events.index = names
        self._bg.index = names

    @property
    def colnames(self):
        return numpy.asarray(self._events.columns)

    @colnames.setter
    def colnames(self, names):
        self._events.columns = names
        self._bg.columns = names

    @property
    def shape(self):
        return self._events.shape

    def __repr__(self):
        return repr(self._events)


def row_stack(matrices):
    """
    Stack DiscoverMatrix objects row-wise.

    Parameters
    ----------
    matrices : sequence of DiscoverMatrix
        Sequence containing the DiscoverMatrix objects to be stacked. The
        matrices must have the same number of columns, and the column names
        must match.

    Returns
    -------
    result : DiscoverMatrix
        The matrix formed by stacking the given matrices.
    """
    assert all((x._events.columns == matrices[0]._events.columns).all() for x in matrices)
    return DiscoverMatrix(
        pandas.concat([x._events for x in matrices], 0),
        pandas.concat([x._bg for x in matrices], 0))
