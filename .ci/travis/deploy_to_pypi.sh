if [[ "$TRAVIS_TAG" == v* ]]; then
  echo "Testing PW=$bob"
  echo "Uploading to PyPI."
  twine upload --skip-existing -u rbarnes -p$PYPI_PASSWORD wheelhouse/*manylinux*whl
  twine upload --skip-existing -u rbarnes -p$PYPI_PASSWORD wheelhouse/*.tar.gz
  echo "Done."
  # if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
  #   echo "Creating a source distribution."
  #   python setup.py sdist
  # else
  #   echo "Creating a binary wheel distribution."
  #   python setup.py bdist_wheel
  # fi
else
  echo "Not deploying."
  echo "Tag is... $TRAVIS_TAG"
  echo "Branch name is... $TRAVIS_BRANCH"
fi