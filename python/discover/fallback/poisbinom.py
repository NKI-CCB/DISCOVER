import numpy


try:
    range = xrange
except NameError:
    # This is Python 3
    pass


def poisbinom(p, K):
    if K > len(p):
        return 0.0, 1.0
    
    logp = numpy.log(p)
    lognotp = numpy.log(1 - p)

    memory = numpy.empty((2, len(p)))

    #k = 0
    memory[0, 0] = lognotp[0]
    for j in xrange(1, len(p)):
        memory[0, j] = lognotp[j] + memory[0, j - 1]

    prev = 0
    curr = 1

    cdf = memory[0, len(p) - 1]

    for k in xrange(1, K + 1):
        if k > 1:
            memory[curr, k - 1] = logp[k - 1] + memory[prev, k - 2]
        else:
            memory[curr, k - 1] = logp[k - 1]
        
        for j in xrange(k, len(p)):
            memory[curr, j] = numpy.logaddexp(
                lognotp[j] + memory[curr, j - 1],
                logp[j] + memory[prev, j - 1])

        cdf = numpy.logaddexp(cdf, memory[curr, len(p) - 1])

        curr, prev = prev, curr

    return numpy.exp(memory[prev, len(p) - 1]), numpy.exp(cdf)


def cdf(p, x):
    return poisbinom(p, x)[1]
