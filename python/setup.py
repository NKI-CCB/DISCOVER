from __future__ import print_function

import glob
import os
import re
import sys

from numpy.distutils import fcompiler
from numpy.distutils import system_info
from numpy.distutils.core import Extension, setup


def scan_regexp(filename, regexp):
    regexp = re.compile(regexp)
    matches = set()

    with open(filename) as stream:
        for line in stream:
            if line[0] not in ["c", "C", "*", "d", "D"]:
                match = regexp.search(line.partition("!")[0])
                if match is not None:
                    matches.add(match.group(1))

    return matches


def topological_sort(graph):
    start_nodes = [node for node in graph if len(graph[node]) == 0]
    graph = { node : graph[node].copy() for node in graph if len(graph[node]) > 0 }
    result = []

    while len(start_nodes) > 0:
        node = start_nodes.pop()
        result.append(node)
        for node2, deps in list(graph.items()):
            if node in deps:
                deps.remove(node)
            if len(deps) == 0:
                del graph[node2]
                start_nodes.append(node2)

    if len(graph) > 0:
        raise ValueError("Input contains circular dependencies")
    else:
        return result


def get_sorted_source_files(src_dir):
    f90_sources = glob.glob(os.path.join(src_dir, "*.f90"))
    f77_sources = glob.glob(os.path.join(src_dir, "*.f"))

    module_sources = { module : filename
                       for filename in f90_sources
                       for module in scan_regexp(filename, r"module (\w+)") }

    dependency_graph = { filename: { module_sources[module]
                                     for module
                                     in scan_regexp(filename, r"use (\w+)") }
                         for filename
                         in f90_sources }
    
    return f77_sources + topological_sort(dependency_graph)


src_dir = "src"
source_files = get_sorted_source_files(src_dir)


compiler = fcompiler.get_default_fcompiler()
for arg in sys.argv:
    if "--fcompiler" in arg:
        compiler = arg[arg.index("=") + 1:]
print("Using %s compiler" % compiler, file=sys.stderr)


if compiler == "gnu95":
    f77flags = []
    f90flags = ["-ffree-line-length-none"]
    opt = ["-Ofast", "-march=native", "-funroll-loops", "-flto", "-fopenmp"]
elif compiler.startswith("intel"):
    f77flags = []
    f90flags = ["-assume realloc_lhs"]
    opt = ["-fast", "-openmp"]
else:
    print("Compiler %s not supported. Let's hope things just work." % compiler, file=sys.stderr)
    f77flags = []
    f90flags = []
    opt = []



blas_info = system_info.get_info("blas_opt")

if len(blas_info) == 0:
    print("No optimised BLAS implementation found. Make sure libopenblas is in the current directory.", file=sys.stderr)
    blas_info = {
        "library_dirs": ["."],
        "libraries": ["openblas"]
    }


if sys.platform.startswith("darwin"):
    link_flags = ["-undefined", "dynamic_lookup", "-bundle"]
else:
    link_flags = []


discover_ext = Extension(
    name="_discover",
    sources=["discover.pyf"] + source_files,
    library_dirs=blas_info["library_dirs"],
    libraries=blas_info["libraries"],
    extra_f90_compile_args=f90flags + opt,
    extra_f77_compile_args=f77flags + opt,
    extra_link_args=opt + link_flags)


if __name__ == "__main__":
    setup(name="discover",
          version="0.9.4",
          description="DISCOVER mutual exclusivity and co-occurrence analysis",
          author="Sander Canisius",
          author_email="s.canisius@nki.nl",
          license="Apache License 2.0",
          ext_modules=[discover_ext],
          packages=["discover",
                    "discover.datasets",
                    "discover.fallback"],
          package_data={"discover.datasets": ["*.bz2"]},
          classifiers=[
              "Intended Audience :: Science/Research",
              "License :: OSI Approved :: Apache Software License",
              "Programming Language :: Python :: 2",
              "Programming Language :: Python :: 3"
          ])
