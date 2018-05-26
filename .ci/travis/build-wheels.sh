#!/bin/bash
set -e -x

mkdir -p /io/wheelhouse

#Ensure the build directory is empty
rm -rf /io/wheelhouse/*

#Make a source distribution
/opt/python/cp36-cp36m/bin/python3.6 /io/wrappers/pyrichdem/setup.py sdist -d /io/wheelhouse/

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r /io/wrappers/pyrichdem/requirements.txt
    #"${PYBIN}/pip" wheel /io/wrappers/pyrichdem/ -w wheelhouse/ #Original many linux command

    #The following is my hackery
    cd /io/wrappers/pyrichdem
    "${PYBIN}/python" setup.py bdist_wheel -d /io/wheelhouse/
    rm -rf build/ richdem.egg-info/ dist/
done

cd /io/wrappers/pyrichdem

# Bundle external shared libraries into the wheels
for whl in /io/wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install richdem --no-index -f /io/wheelhouse
    #(cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
done
