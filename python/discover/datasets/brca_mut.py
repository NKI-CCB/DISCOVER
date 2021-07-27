import bz2
import numpy
import pandas
import os


DATA_FILE_NAME = os.path.join(os.path.dirname(__file__), "brca_mut.npy.bz2")


def load():
    with bz2.BZ2File(DATA_FILE_NAME) as stream:
        values = numpy.load(stream)
        index = numpy.load(stream).astype(str)
        columns = numpy.load(stream).astype(str)

    return pandas.DataFrame(values, index=index, columns=columns)
