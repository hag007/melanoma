!#/bin/bash
DIRECTORY = ../../../python27_test

if [ ! -d "$DIRECTORY" ]; then
  mkdir $DIRECTORY
fi

virtualenv $DIRECTORY
source $DIRECTORY/bin/activate.csh
pip install -r dependencies.txt
deactivate