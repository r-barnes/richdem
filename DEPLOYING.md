1. Update CHANGELOG.md
2. Bump version in `include/richdem/common/version.hpp`
3. Bump version in `wrappers/pyrichdem/setup.py`


Setting up a Docker environment
===============================

Acquire docker:

    sudo apt install docker.io

Give user docker access:

    sudo usermod -a -G docker $USER

Login and logout of user for permissions to take effect.

Acquire docker images:

    docker pull quay.io/pypa/manylinux1_x86_64
    docker pull quay.io/pypa/manylinux1_i686

List docker images:

    docker images

Run docker images:

    docker run -i -t 154333521e84 /bin/bash

Get RichDEM:

    git clone https://github.com/r-barnes/richdem.git
    mv richdem io

Compile compile compile:

    cd /io/travis
    ./build-wheels.sh

Copy built files from docker:

    docker cp <containerId>:/io/wheelhouse/* /z/

    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp27-cp27m-linux_i686.whl               /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp33-cp33m-linux_i686.whl              /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp35-cp35m-linux_i686.whl       /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp27-cp27m-manylinux1_i686.whl          /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp33-cp33m-manylinux1_i686.whl         /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp35-cp35m-manylinux1_i686.whl       /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp27-cp27mu-linux_i686.whl              /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp34-cp34m-linux_i686.whl              /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp36-cp36m-linux_i686.whl       /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp27-cp27mu-manylinux1_i686.whl         /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp34-cp34m-manylinux1_i686.whl         /z/
    docker cp 23659cfa50b3:/io/wheelhouse/richdem-0.2.0-cp36-cp36m-manylinux1_i686.whl       /z/

    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp27-cp27m-linux_x86_64.whl         /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp34-cp34m-linux_x86_64.whl /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp27-cp27m-manylinux1_x86_64.whl    /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp34-cp34m-manylinux1_x86_64.whl /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp27-cp27mu-linux_x86_64.whl        /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp35-cp35m-linux_x86_64.whl /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp27-cp27mu-manylinux1_x86_64.whl   /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp35-cp35m-manylinux1_x86_64.whl /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp33-cp33m-linux_x86_64.whl         /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp36-cp36m-linux_x86_64.whl /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp33-cp33m-manylinux1_x86_64.whl    /z/
    docker cp f7f25ab807da:/io/wheelhouse/richdem-0.2.0-cp36-cp36m-manylinux1_x86_64.whl /z/

List docker containers:

    docker ps -a

Remove a stopped container:

    docker rm <CONTAINER ID>

Remove a docker image:

    docker rmi <IMAGE ID>

Start a stopped container:

    docker start -i <CONTAINER ID>



Compiling a source distribution
===============================

    python3 setup.py sdist



Uploading packages to PyPI
==========================

    twine upload dist/richdem-0.3.0.tar.gz
    twine upload /z/richdem-0.2.0-cp*-manylinux*



Create a Zenodo Release
==========================

Go [here](https://github.com/r-barnes/richdem/releases/new) and fill out the form for the new tag.