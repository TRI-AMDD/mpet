#!/bin/bash
#Add test dir to environment
export PYTHONPATH=$PYTHONPATH:$( dirname "${BASH_SOURCE[0]}" )/../tests

#pass through to test script
python ./bin/run_tests.py --test_dir $( dirname "${BASH_SOURCE[0]}" )/../tests $@
