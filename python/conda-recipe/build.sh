if [[ $(uname) == 'Darwin' ]]; then
    export LDFLAGS="$LDFLAGS -undefined dynamic_lookup"
else
    export LDFLAGS="$LDFLAGS -shared"
fi

export F90=$GFORTRAN

$PYTHON setup.py install
