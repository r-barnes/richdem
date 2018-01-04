Release Steps
==================================================

Updating Documentation
--------------------------------------------------

1. Enter the `docs` directory.

2. Run `make html`

3. Enter the `docs/richdem-docs` directory. Commit changes.

4. Go back to the `docs` directory. Commit changes.

5. Push.

6. This will trigger a ReadTheDocs build and a Travis build.


Updating Wheels
--------------------------------------------------

Acquire the manylinux build image with:

    docker pull quay.io/pypa/manylinux1_x86_64

Start the docker manylinux build image with:

    docker run -i -t e0e55980c200 /bin/bash

You can list running images with:

    docker ps

From within the container run the `travis/build-wheels.sh` script to generate
documentation.

From without the container run, e.g.:

    docker cp mystifying_darwin:/io/wheelhouse/richdem-0.0.9-cp27-cp27m-linux_x86_64.whl ./

Upload the wheel to PyPI with

    twine upload richdem-0.0.9-cp*


Updating Source Dist
--------------------------------------------------

To generate a source distribution, run

    python3 setup.py sdist

Upload it to PyPI with, e.g.

    twine upload dist/richdem-0.0.10.tar.gz