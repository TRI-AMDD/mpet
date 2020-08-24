#!/bin/bash
set -x
set -e
python mpettest.py | tee test_output.txt
set +e
grep fail test_output.txt
if [ "$?" -eq 0 ]; then 
  echo "Testing Failed"
  exit 1
else
  echo "All Tests Passed"
  exit 0
fi
