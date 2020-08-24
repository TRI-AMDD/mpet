#!/bin/bash
set -e

python mpettest.py | tee test_output.txt
grep fail test_output.txt > /dev/null
if [ "$?$" -eq 1 ] then 
  echo "Testing Failed"
  exit 1
else
  echo "All Tests Passed"
  exit 0
fi
