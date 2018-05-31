# if [[ "$TRAVIS_TAG" == v* ]]; then
  echo "Testing PW=$bob"
  echo "Uploading to PyPI."
  twine upload -u rbarnes -p$PYPI_PASSWORD wheelhouse/*
  echo "Done."
  # if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
  #   echo "Creating a source distribution."
  #   python setup.py sdist
  # else
  #   echo "Creating a binary wheel distribution."
  #   python setup.py bdist_wheel
  # fi
# else
#   echo "Not deploying."
#   echo "Tag is... $TRAVIS_TAG"
#   echo "Branch name is... $TRAVIS_BRANCH"
# fi