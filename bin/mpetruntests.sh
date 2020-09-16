#!/bin/bash
set -x
set -e
rm -rf workdir
mkdir workdir
cd workdir

cp ../run_tests.py .
ln -s ../../mpet .
ln -s ../../tests .

#run tests for current commit
coverage run --source=../../mpet/ run_tests.py ../../bin/workdir/modified > /dev/null
COVERALLS_REPO_TOKEN=$COVERALLS_REPO_TOKEN coveralls || true #Dont fret if if fails

#run tests for stable branch
git fetch origin $1:$1 #get stable branch before cloning
git clone ../../ stable_branch/
cd stable_branch/
git checkout $1
cd ../
ln -sf stable_branch/mpet .
python run_tests.py ../../bin/workdir/stable > /dev/null

