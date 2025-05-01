project(
  'discover',
  ['c', 'fortran'],
  version: '0.9.6',
  license: 'Apache-2.0',
  default_options: ['b_lto=true']
)


# Set specific compiler options

fc = meson.get_compiler('fortran')

arch = get_option('arch')
if arch != ''
  add_project_arguments(fc.get_supported_arguments('-march=' + arch), language: 'fortran')
endif

if fc.get_id() == 'gcc'
  add_project_arguments('-ffree-line-length-none', language: 'fortran')
  add_project_arguments('-ftree-vectorize', '-funroll-loops', '-finline-functions', language: 'fortran')
  
  if get_option('fast')
    add_project_arguments('-Ofast', language: 'fortran')

    if not get_option('stack-arrays')
      add_project_arguments('-fno-stack-arrays', language: 'fortran')
    endif
  elif get_option('stack-arrays')
    add_project_arguments('-fstack-arrays', language: 'fortran')
  endif
elif fc.get_id() == 'intel'
  add_project_arguments('-assume realloc_lhs', language: 'fortran')

  if get_option('fast')
    add_project_arguments('-fast', language: 'fortran')
  endif

  # see heap-arrays https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/heap-arrays.html
endif


py_mod = import('python')
py = py_mod.find_installation(pure: false)


# Find numpy and f2py include directories

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)


# Find BLAS library dependencies

blas_name = get_option('blas')

if blas_name == 'auto'
  numpy_version = run_command(py, ['-c', 'import numpy; print(numpy.__version__)'], check: true).stdout().strip()

  if numpy_version.version_compare('>=1.25.0')
    blas_name = run_command(py,
			    ['-c', 'import numpy; print(numpy.show_config("dicts")["Build Dependencies"]["blas"]["name"])'],
			    check: true
			   ).stdout().strip()
    blas_deps = dependency(blas_name)
  else
    blas_libs = run_command(py,
			    ['-c', 'import numpy; print(*numpy.__config__.get_info("blas_opt")["libraries"], sep="\\n")'],
			    check : true
			   ).stdout().splitlines()

    blas_dirs = run_command(py,
			    ['-c', 'import numpy; print(*numpy.__config__.get_info("blas_opt")["library_dirs"], sep="\\n")'],
			    check : true
			   ).stdout().splitlines()

    blas_deps = []
    foreach lib : blas_libs
      blas_deps += fc.find_library(lib, dirs: blas_dirs)
    endforeach
  endif
else
  blas_deps = dependency(blas_name)
endif


# Generate f2py wrappers

f2py = find_program('f2py')
f2py_wrappers = custom_target(
  'f2py-wrappers',
  input: 'discover.pyf',
  output: ['_discover-f2pywrappers2.f90', '_discovermodule.c'],
  command: [f2py, '@INPUT@']
)


# Build the extension module

py_dep = py.dependency()
openmp_dep = dependency('openmp', required: get_option('openmp'))

py.extension_module(
  '_discover',
  [
    'src/types.f90',
    'src/logspace.f90',
    'src/util.f90',
    'src/uniinv.f90',
    'src/fdr.f90',
    'src/lbfgsb.f',
    'src/linpack.f',
    'src/maxent.f90',
    'src/poisbinom.f90',
    'src/timer.f',
  ],
  f2py_wrappers,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : [py_dep, blas_deps, openmp_dep],
  install : true
)


# Copy Python files

py.install_sources(
  [
    'discover/__init__.py',
    'discover/background.py',
    'discover/data.py',
    'discover/grouptest.py',
    'discover/pairwise.py',
    'discover/plots.py',
    'discover/stats.py',
    'discover/util.py',
    'discover/datasets/__init__.py',
    'discover/datasets/brca_mut.py',
    'discover/datasets/brca_mut.npy.bz2',
    'discover/fallback/__init__.py',
    'discover/fallback/poisbinom.py'
  ],
  preserve_path: true
)
