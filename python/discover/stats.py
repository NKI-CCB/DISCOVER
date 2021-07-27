import numpy


def false_discovery_rate(p, pi0=1.0):
    if not 0 <= pi0 <= 1:
        raise ValueError("Invalid value for pi0: %s. Legal values are between 0 and 1" % pi0)

    nna = ~numpy.isnan(p)
    q = numpy.full_like(p, numpy.nan)
    p = p[nna]

    i = numpy.arange(len(p), 0, -1)
    o = p.argsort()[::-1]
    ro = o.argsort()
    q[nna] = numpy.minimum(1, numpy.minimum.accumulate(float(pi0) * len(p) / i * p[o])[ro])
    
    return q
