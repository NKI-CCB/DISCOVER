import importlib
import pkgutil


def get_dataset_names():
    return [name for importer, name, ispkg in pkgutil.iter_modules(__path__)]


def load_dataset(name):
    dataset_module = importlib.import_module(".%s" % name, __package__)
    return dataset_module.load()
