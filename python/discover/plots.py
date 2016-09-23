import matplotlib
import matplotlib.pyplot as plt
import numpy


def eventPlot(events, tickDist=None):
    geneOrder = numpy.asarray(events).sum(1).argsort()
    sampleOrder = numpy.lexsort(numpy.asarray(events)[geneOrder])[::-1]

    plt.imshow(
        numpy.asarray(events)[geneOrder[:, numpy.newaxis], sampleOrder],
        aspect="auto", interpolation="none", cmap=matplotlib.colors.ListedColormap(["#dddddd", "#377eb8"]),
        origin="lower")

    if tickDist is not None:
        plt.xticks(numpy.arange(0, events.shape[1], tickDist), [])
    
    plt.grid(ls="-", c="white", axis="x")
    plt.setp(list(plt.gca().spines.values()), color="#888888", alpha=0)

    for tic in plt.gca().xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    for tic in plt.gca().yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    try:
        # If `events` is a pandas DataFrame, we use its index as row labels
        plt.yticks(numpy.arange(len(events)), events.index[geneOrder])
    except AttributeError:
        pass

    plt.gca().yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(
        numpy.linspace(-0.5, len(events) - 0.5, len(events) + 1)))
    plt.grid(ls="-", c="white", axis="y", lw=4, alpha=1, which="minor")
    plt.grid(ls="None", axis="y", which="major")
    for tic in plt.gca().yaxis.get_minor_ticks():
        tic.tick1On = tic.tick2On = False

    plt.gcf().patch.set_facecolor("white")
