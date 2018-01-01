#!/bin/bash
set -e -x

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r /io/wrappers/pyrichdem/requirements.txt
    "${PYBIN}/pip" wheel /io/wrappers/pyrichdem/ -w wheelhouse/

    cd /io/wrappers/pyrichdem
    "${PYBIN}/python" setup.py bdist_wheel -d ../../wheelhouse/
    rm -rf build/ richdem.egg-info/ dist/
done

cd /io/

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install richdem --no-index -f /io/wheelhouse
    #(cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
done
